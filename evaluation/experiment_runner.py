import os
import copy
import numpy as np
import torch
from evaluation.metrics import compute_all_metrics
from evaluation.logger import ExperimentLogger
from evaluation.plotting import plot_optimization_curve, plot_all

NUM_RUNS = 3


def decode_to_smiles(token_indices, vocab):
    idx_to_char = {v: k for k, v in vocab.items()}
    tokens = token_indices.squeeze(0).tolist()
    smiles = ""
    for t in tokens:
        ch = idx_to_char.get(t, "")
        if ch in ("<eos>", "<pad>", ""):
            break
        smiles += ch
    return smiles


def generate_with_gradient(vae, predictor, objective_fn, device, n_samples=100, n_steps=50, lr=0.05):
    vae.eval()
    predictor.eval()

    z_list = []
    trajectory = []

    for _ in range(n_samples):
        z = torch.randn(1, vae.latent_dim, device=device, requires_grad=True)
        optimizer = torch.optim.Adam([z], lr=lr)

        best_score = -float("inf")
        for step in range(n_steps):
            optimizer.zero_grad()
            score = objective_fn(z, predictor)
            loss = -score
            loss.backward()
            optimizer.step()
            best_score = max(best_score, score.item())
            if _ == 0:
                trajectory.append(score.item())

        z_list.append(z.detach())

    return torch.cat(z_list, dim=0), trajectory


def generate_with_bayesian(vae, dataloader, predictor, objective_fn, device):
    from bayesian_optimization import run_bayesian_optimization

    molecule, best_z, result = run_bayesian_optimization(
        vae=vae,
        dataloader=dataloader,
        predictor=predictor,
        objective_fn=objective_fn,
        device=device,
    )
    trajectory = [-v for v in result.func_vals.tolist()]
    best_z_tensor = torch.tensor(best_z, dtype=torch.float).to(device).view(1, -1)
    return best_z_tensor, trajectory


def run_single(vae_instance, opt_name, predictor, objective_fn, dataloader, vocab,
               training_smiles_set, device, n_samples=100):
    if opt_name == "gradient":
        z_all, trajectory = generate_with_gradient(
            vae_instance, predictor, objective_fn, device, n_samples=n_samples
        )
    else:
        z_best, trajectory = generate_with_bayesian(
            vae_instance, dataloader, predictor, objective_fn, device
        )
        z_all = z_best.expand(n_samples, -1)

    smiles_list = []
    with torch.no_grad():
        for i in range(z_all.shape[0]):
            z_i = z_all[i].unsqueeze(0)
            token_indices = vae_instance.decode(z_i, device)
            smi = decode_to_smiles(token_indices, vocab)
            if smi:
                smiles_list.append(smi)

    metrics = compute_all_metrics(
        smiles_list=smiles_list,
        training_smiles_set=training_smiles_set,
        predictor=predictor,
        device=device,
    )
    metrics["trajectory"] = trajectory
    return metrics


def run_experiments(
    models,
    vae,
    predictor,
    objective_fn,
    dataloader,
    vocab,
    training_smiles_set,
    device,
    n_samples=100,
    k=NUM_RUNS,
    save_dir="evaluation",
):
    """
    models: dict[model_name] -> encoder module (or None to use vae.encoder as-is)

    Runs each (model, optimizer) config k times.
    Aggregates mean and std across runs.
    """
    logger = ExperimentLogger()
    optimizer_names = ["gradient", "bayesian"]

    for model_name, encoder in models.items():
        for opt_name in optimizer_names:
            for run_id in range(k):
                print(f"\n[RUN] model={model_name} | optimizer={opt_name} | run={run_id + 1}/{k}")

                # Model isolation: deep copy VAE and swap encoder
                vae_instance = copy.deepcopy(vae)
                if encoder is not None:
                    vae_instance.encoder = copy.deepcopy(encoder)

                try:
                    metrics = run_single(
                        vae_instance=vae_instance,
                        opt_name=opt_name,
                        predictor=predictor,
                        objective_fn=objective_fn,
                        dataloader=dataloader,
                        vocab=vocab,
                        training_smiles_set=training_smiles_set,
                        device=device,
                        n_samples=n_samples,
                    )

                    config = {"model": model_name, "optimizer": opt_name}
                    logger.log_run(config, metrics, run_id=run_id)

                    plot_optimization_curve(
                        metrics.get("trajectory", []),
                        model_name,
                        f"{opt_name}_run{run_id}",
                        save_dir=os.path.join(save_dir, "plots"),
                    )

                except Exception as e:
                    print(f"[ERROR] {model_name}/{opt_name}/run{run_id} failed: {e}")

    logger.aggregate()
    log_path = os.path.join(save_dir, "logs", "results.json")
    logger.save(log_path)

    store = logger.get_store()
    plot_all(store, save_dir=os.path.join(save_dir, "plots"))

    return store
