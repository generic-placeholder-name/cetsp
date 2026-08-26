#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_manifest(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty snapshot manifest: {path}")
    return rows


def parse_run(value):
    if "=" not in value:
        raise argparse.ArgumentTypeError("run must be LABEL=MANIFEST.csv")
    label, path = value.split("=", 1)
    if not label or not path:
        raise argparse.ArgumentTypeError("run must be LABEL=MANIFEST.csv")
    return label, Path(path)


def series(rows, column):
    return [float(row[column]) for row in rows]


def plot(runs, output, title):
    figure, axes = plt.subplots(
        3,
        1,
        figsize=(12, 11),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [1.25, 1.0, 0.55]},
    )
    colors = plt.colormaps["tab10"]

    for index, (label, path) in enumerate(runs):
        rows = read_manifest(path)
        color = colors(index)
        generations = series(rows, "generation")

        axes[0].plot(
            generations,
            series(rows, "best_cost"),
            color=color,
            linewidth=2.2,
            label=f"{label} best",
        )
        axes[0].plot(
            generations,
            series(rows, "mean_cost"),
            color=color,
            linewidth=1.5,
            linestyle="--",
            label=f"{label} mean",
        )
        axes[0].plot(
            generations,
            series(rows, "worst_cost"),
            color=color,
            linewidth=1.0,
            linestyle=":",
            alpha=0.8,
            label=f"{label} worst",
        )

        axes[1].plot(
            generations,
            series(rows, "mean_pairwise_edit_distance"),
            color=color,
            linewidth=2.2,
            label=f"{label} mean pairwise",
        )
        axes[1].plot(
            generations,
            series(rows, "min_pairwise_edit_distance"),
            color=color,
            linewidth=1.5,
            linestyle="--",
            label=f"{label} minimum pairwise",
        )

        axes[2].step(
            generations,
            series(rows, "population"),
            color=color,
            linewidth=1.8,
            where="post",
            label=label,
        )

    axes[0].set_ylabel("Tour length")
    axes[0].set_title("Population quality")
    axes[0].legend(ncols=3, fontsize=8)

    axes[1].set_ylabel("Edit distance (%)")
    axes[1].set_title("Ordering diversity")
    axes[1].legend(ncols=2, fontsize=8)

    axes[2].set_ylabel("Members")
    axes[2].set_xlabel("Generation")
    axes[2].set_title("Population expansion and culling")
    axes[2].legend(fontsize=8)

    for axis in axes:
        axis.grid(color="#d1d5db", linewidth=0.7, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)

    figure.suptitle(title, fontsize=15)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description="Plot MA-CETSP population snapshot histories")
    parser.add_argument("output", type=Path)
    parser.add_argument("runs", nargs="+", type=parse_run)
    parser.add_argument("--title", default="MA-CETSP population history")
    args = parser.parse_args()
    plot(args.runs, args.output, args.title)
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
