from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base_dir = Path(__file__).resolve().parents[1]

summary_path = base_dir / "results" / "tables" / "pbmc10k_decontx_quick_benchmark.csv"
figure_dir = base_dir / "results" / "figures"
figure_dir.mkdir(parents=True, exist_ok=True)

summary = pd.read_csv(summary_path)

df = summary[summary["method"] == "DecontX"].copy()

def condition_order_key(condition):
    condition = str(condition)

    if condition.startswith("1gene"):
        ambient_rank = 0
    elif condition.startswith("100genes"):
        ambient_rank = 1
    else:
        ambient_rank = 2

    if "5pct" in condition and "15pct" not in condition:
        contam_rank = 0
    elif "15pct" in condition:
        contam_rank = 1
    elif "30pct" in condition:
        contam_rank = 2
    else:
        contam_rank = 3

    return ambient_rank, contam_rank

df["order"] = df["condition"].apply(condition_order_key)
df = df.sort_values("order")

x = np.arange(len(df))
width = 0.25

plt.figure(figsize=(10, 5))

plt.bar(
    x - width,
    df["designed_contamination_fraction"],
    width,
    label="Designed contamination"
)

plt.bar(
    x,
    df["mean_decontx_estimated_contamination"],
    width,
    label="Mean DecontX estimate"
)

plt.bar(
    x + width,
    df["median_decontx_estimated_contamination"],
    width,
    label="Median DecontX estimate"
)

plt.xticks(x, df["condition"], rotation=45, ha="right")
plt.ylabel("Contamination fraction")
plt.xlabel("Condition")
plt.title("Designed contamination vs DecontX estimated contamination")
plt.legend()
plt.tight_layout()

out_path = figure_dir / "pbmc10k_decontx_estimated_vs_designed_contamination.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight")
plt.close()

print("Saved figure:")
print(out_path)
