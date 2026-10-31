"""Draw the README figures from benchmarks/results/*.json."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "benchmarks" / "results"
FIGS = ROOT / "figures"


def main() -> None:
    FIGS.mkdir(exist_ok=True)
    k15 = json.loads((RESULTS / "mapper_benchmark.json").read_text())
    k11_path = RESULTS / "mapper_benchmark_k11.json"

    fig, ax = plt.subplots(figsize=(6, 4))
    for data, label in [(k15, "k = 15")] + (
        [(json.loads(k11_path.read_text()), "k = 11")] if k11_path.exists() else []
    ):
        rates = [100 * r["error_rate"] for r in data["mapper"]]
        ax.plot(rates, [r["correct_pct"] for r in data["mapper"]], marker="o", label=label)
    ax.set_xlabel("Simulated per-base error rate (%)")
    ax.set_ylabel("Reads mapped to the correct position (%)")
    ax.set_title(f"Mapper accuracy on SARS-CoV-2 ({k15['reference_bp']:,} bp), 150 bp reads")
    ax.set_ylim(80, 101)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGS / "mapper_accuracy.png", dpi=150)

    dp = k15["dp_scaling"]
    fig, ax = plt.subplots(figsize=(6, 4))
    cells = [r["cells"] for r in dp]
    secs = [r["seconds"] for r in dp]
    ax.loglog(cells, secs, marker="o", label="measured (affine global, pure Python)")
    slope = secs[-1] / cells[-1]
    ax.loglog(cells, [c * slope for c in cells], "--", label="linear in cells = O(n*m)")
    ax.set_xlabel("DP table cells (n+1)*(m+1)")
    ax.set_ylabel("Seconds")
    ax.set_title("Needleman-Wunsch/Gotoh runtime grows with n*m")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGS / "dp_scaling.png", dpi=150)
    print("wrote", sorted(p.name for p in FIGS.iterdir()))


if __name__ == "__main__":
    main()
