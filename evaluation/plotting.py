import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_optimization_curve(objective_values, model_name, optimizer_name, save_dir="evaluation/plots"):
    os.makedirs(save_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(range(1, len(objective_values) + 1), objective_values, marker="o", linewidth=2, markersize=4)
    ax.set_xlabel("Iteration", fontsize=13)
    ax.set_ylabel("Objective Value", fontsize=13)
    ax.set_title(f"Optimization Curve — {model_name} / {optimizer_name}", fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    fname = os.path.join(save_dir, f"opt_curve_{model_name}_{optimizer_name}.png")
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"[PLOT] Saved optimization curve → {fname}")


def plot_model_comparison(store, metric="mean_qed", save_dir="evaluation/plots"):
    """
    store: dict[model][optimizer] = {"mean": {...}, "std": {...}, "runs": [...]}
    Grouped bar chart with error bars.
    """
    os.makedirs(save_dir, exist_ok=True)

    model_names = list(store.keys())
    optimizer_names = list(next(iter(store.values())).keys())

    x = list(range(len(model_names)))
    n_opts = len(optimizer_names)
    bar_width = 0.35
    fig, ax = plt.subplots(figsize=(9, 6))

    for i, opt in enumerate(optimizer_names):
        mean_values = [
            store[m].get(opt, {}).get("mean", {}).get(metric, 0.0)
            for m in model_names
        ]
        std_values = [
            store[m].get(opt, {}).get("std", {}).get(metric, 0.0)
            for m in model_names
        ]
        offsets = [xi + (i - n_opts / 2 + 0.5) * bar_width for xi in x]
        bars = ax.bar(
            offsets,
            mean_values,
            width=bar_width,
            yerr=std_values,
            capsize=5,
            label=opt,
            error_kw={"elinewidth": 1.5, "ecolor": "black"},
        )
        for bar, val in zip(bars, mean_values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.02,
                f"{val:.3f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(model_names, fontsize=12)
    ax.set_ylabel(metric.replace("_", " ").title(), fontsize=13)
    ax.set_title(f"Model Comparison — {metric.replace('_', ' ').title()}", fontsize=14)
    ax.legend(fontsize=11)
    ax.set_ylim(0, 1.15)
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()

    fname = os.path.join(save_dir, f"comparison_{metric}.png")
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"[PLOT] Saved comparison bar chart → {fname}")


def plot_all(store, save_dir="evaluation/plots"):
    plot_model_comparison(store, metric="mean_qed", save_dir=save_dir)
    plot_model_comparison(store, metric="success_rate", save_dir=save_dir)
