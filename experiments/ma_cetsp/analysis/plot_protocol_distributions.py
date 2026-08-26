#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_rows(path):
    with path.open(newline="", encoding="utf-8-sig") as input_file:
        return list(csv.DictReader(input_file))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--runs",
        type=Path,
        default=Path("results/mennell_hybrid/tables/mennell_protocol_runs.csv"),
    )
    parser.add_argument(
        "--lei-hao",
        type=Path,
        default=Path("results/mennell_hybrid/tables/lei_hao_run_statistics.csv"),
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("results/mennell_hybrid/tables/mennell_protocol_summary.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("paper/img/bubbles_hybrid_distribution.png"),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    names = [f"bubbles{index}.txt" for index in range(4, 10)]
    positions = np.arange(len(names), dtype=float)

    references = {row["instance"]: row for row in read_rows(args.lei_hao)}
    summaries = {row["instance"]: row for row in read_rows(args.summary)}
    values = {name: [] for name in names}
    for row in read_rows(args.runs):
        if row["instance"] in values:
            mean = float(references[row["instance"]]["lei_hao_mean"])
            values[row["instance"]].append(100.0 * (float(row["final"]) / mean - 1.0))

    figure, axis = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)
    boxes = axis.boxplot(
        [values[name] for name in names],
        positions=positions,
        widths=0.52,
        patch_artist=True,
        medianprops={"color": "#174a7e", "linewidth": 1.8},
        boxprops={"facecolor": "#8ecae6", "edgecolor": "#174a7e"},
        whiskerprops={"color": "#174a7e"},
        capprops={"color": "#174a7e"},
        flierprops={"marker": ""},
    )
    boxes["boxes"][0].set_label("Hybrid distribution")

    rng = np.random.default_rng(20260825)
    for position, name in zip(positions, names):
        jitter = rng.uniform(-0.16, 0.16, len(values[name]))
        axis.scatter(
            position + jitter,
            values[name],
            color="#174a7e",
            alpha=0.58,
            s=15,
            linewidth=0,
        )

    lei_sd = []
    lei_best = []
    for name in names:
        mean = float(references[name]["lei_hao_mean"])
        lei_sd.append(100.0 * float(references[name]["lei_hao_sd"]) / mean)
        lei_best.append(100.0 * (float(summaries[name]["lei_hao"]) / mean - 1.0))

    axis.errorbar(
        positions + 0.27,
        np.zeros(len(names)),
        yerr=lei_sd,
        fmt="o",
        color="#d1495b",
        capsize=4,
        markersize=4,
        label=r"Lei--Hao mean $\pm$ SD",
    )
    axis.scatter(
        positions + 0.27,
        lei_best,
        marker="D",
        color="#6a040f",
        s=28,
        label="Lei--Hao best",
        zorder=4,
    )
    axis.axhline(0.0, color="#666666", linewidth=0.8, linestyle="--")
    axis.set_xticks(positions, [f"B{index}" for index in range(4, 10)])
    axis.set_ylabel("Gap from Lei--Hao mean (%)")
    axis.set_xlabel("Bubbles instance")
    axis.legend(frameon=False, ncols=3, loc="upper left")
    axis.grid(axis="y", color="#dddddd", linewidth=0.7)
    axis.set_axisbelow(True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=220)


if __name__ == "__main__":
    main()
