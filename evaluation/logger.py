import json
import os
from datetime import datetime


class ExperimentLogger:
    def __init__(self):
        # Nested: runs[model][optimizer] = {"runs": [], "mean": {}, "std": {}}
        self.store = {}

    def log_run(self, config, metrics, run_id=0):
        model = config.get("model")
        optimizer = config.get("optimizer")

        if model not in self.store:
            self.store[model] = {}
        if optimizer not in self.store[model]:
            self.store[model][optimizer] = {"runs": []}

        entry = {
            "run_id": run_id,
            "timestamp": datetime.utcnow().isoformat(),
            "config": config,
            "metrics": metrics,
        }
        self.store[model][optimizer]["runs"].append(entry)

        print(
            f"[LOG] run={run_id} | model={model} | optim={optimizer} | "
            f"validity={metrics.get('validity', 0):.3f} | "
            f"qed={metrics.get('mean_qed', 0):.3f} | "
            f"success={metrics.get('success_rate', 0):.3f} | "
            f"diversity={metrics.get('diversity', 0):.3f}"
        )

    def aggregate(self):
        import numpy as np

        for model, optimizers in self.store.items():
            for optimizer, data in optimizers.items():
                runs = data["runs"]
                if not runs:
                    continue

                metric_keys = runs[0]["metrics"].keys()
                mean_metrics = {}
                std_metrics = {}

                for key in metric_keys:
                    values = [r["metrics"].get(key, 0.0) for r in runs]
                    mean_metrics[key] = float(np.mean(values))
                    std_metrics[key] = float(np.std(values))

                self.store[model][optimizer]["mean"] = mean_metrics
                self.store[model][optimizer]["std"] = std_metrics

    def save(self, filepath):
        self.aggregate()
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, "w") as f:
            json.dump(self.store, f, indent=2)
        print(f"[LOG] Results saved to {filepath}")

    def get_store(self):
        return self.store
