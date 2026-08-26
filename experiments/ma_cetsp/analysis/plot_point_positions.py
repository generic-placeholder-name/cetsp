#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
from statistics import mean, median

import matplotlib.pyplot as plt

from plot_seed_tours import read_instance, read_seeds


def parse_run(value):
    if "=" not in value:
        raise argparse.ArgumentTypeError("run must be LABEL=POPULATION.seeds")
    label, path = value.split("=", 1)
    if not label or not path:
        raise argparse.ArgumentTypeError("run must be LABEL=POPULATION.seeds")
    return label, Path(path)


def empirical_cdf(values):
    ordered = sorted(values)
    return ordered, [(index + 1) / len(ordered) for index in range(len(ordered))]


def summarize(circles, tours):
    min_x = min(x - radius for x, _, radius in circles)
    max_x = max(x + radius for x, _, radius in circles)
    min_y = min(y - radius for _, y, radius in circles)
    max_y = max(y + radius for _, y, radius in circles)
    tolerance = 1e-10 * max(1.0, max_x - min_x, max_y - min_y)

    radial_positions = []
    edge_lengths = []
    unique_fractions = []
    zero_edge_fractions = []
    straight_fractions = []

    for tour in tours:
        points = [(x, y) for _, x, y in tour]
        coordinate_keys = {
            (round(x / tolerance), round(y / tolerance)) for x, y in points
        }
        unique_fractions.append(len(coordinate_keys) / len(points))

        tour_edges = []
        for first, second in zip(points, points[1:] + points[:1]):
            tour_edges.append(math.dist(first, second))
        edge_lengths.extend(tour_edges)
        zero_edge_fractions.append(
            sum(length <= tolerance for length in tour_edges) / len(tour_edges)
        )

        straight = 0
        eligible = 0
        for previous, current, following in zip(
            points[-1:] + points[:-1], points, points[1:] + points[:1]
        ):
            incoming = (current[0] - previous[0], current[1] - previous[1])
            outgoing = (following[0] - current[0], following[1] - current[1])
            incoming_length = math.hypot(*incoming)
            outgoing_length = math.hypot(*outgoing)
            if incoming_length <= tolerance or outgoing_length <= tolerance:
                continue
            eligible += 1
            sine = abs(
                incoming[0] * outgoing[1] - incoming[1] * outgoing[0]
            ) / (incoming_length * outgoing_length)
            cosine = (
                incoming[0] * outgoing[0] + incoming[1] * outgoing[1]
            ) / (incoming_length * outgoing_length)
            if sine <= 1e-7 and cosine > 0.0:
                straight += 1
        straight_fractions.append(straight / eligible if eligible else 0.0)

        for node_id, x, y in tour:
            center_x, center_y, radius = circles[node_id]
            if radius > tolerance:
                radial_positions.append(
                    math.hypot(x - center_x, y - center_y) / radius
                )

    positive_edges = [length for length in edge_lengths if length > tolerance]
    return {
        "radial": radial_positions,
        "edges": positive_edges,
        "unique": mean(unique_fractions),
        "zero_edges": mean(zero_edge_fractions),
        "straight": mean(straight_fractions),
        "center": sum(value <= 0.01 for value in radial_positions)
        / len(radial_positions),
        "boundary": sum(value >= 0.99 for value in radial_positions)
        / len(radial_positions),
        "median_edge": median(positive_edges),
    }


def plot(circles, runs, output, title):
    summaries = []
    for label, path in runs:
        node_count, tours = read_seeds(path)
        if node_count != len(circles):
            raise ValueError(
                f"{path}: {node_count} nodes, expected {len(circles)}"
            )
        summaries.append((label, summarize(circles, tours)))

    figure, axes = plt.subplots(1, 3, figsize=(16, 5.2), constrained_layout=True)
    colors = plt.colormaps["tab10"]
    for index, (label, summary) in enumerate(summaries):
        color = colors(index)
        x, y = empirical_cdf(summary["radial"])
        axes[0].plot(x, y, label=label, color=color, linewidth=2.0)
        x, y = empirical_cdf(summary["edges"])
        axes[1].plot(x, y, label=label, color=color, linewidth=2.0)

    axes[0].set_xlabel("Distance from center / circle radius")
    axes[0].set_ylabel("Fraction of nodes")
    axes[0].set_xlim(0.0, 1.01)
    axes[0].set_title("Position within assigned circle")
    axes[0].legend(fontsize=8)

    axes[1].set_xlabel("Distance to next tour point")
    axes[1].set_ylabel("Fraction of nonzero edges")
    axes[1].set_xscale("log")
    axes[1].set_title("Spacing between consecutive points")
    axes[1].legend(fontsize=8)

    names = [label for label, _ in summaries]
    metrics = [
        ("Duplicate\npoints", "unique", True),
        ("Zero-length\nedges", "zero_edges", False),
        ("Straight-through\npoints", "straight", False),
        ("At circle\ncenter", "center", False),
        ("Near circle\nboundary", "boundary", False),
    ]
    width = 0.8 / len(summaries)
    positions = range(len(metrics))
    for run_index, (label, summary) in enumerate(summaries):
        values = [
            (1.0 - summary[key]) if invert else summary[key]
            for _, key, invert in metrics
        ]
        offsets = [
            position - 0.4 + width / 2 + run_index * width
            for position in positions
        ]
        axes[2].bar(
            offsets,
            [100.0 * value for value in values],
            width=width,
            label=label,
            color=colors(run_index),
        )
    axes[2].set_xticks(list(positions), [label for label, _, _ in metrics])
    axes[2].set_ylabel("Nodes or edges (%)")
    axes[2].set_title("Positional structure")
    axes[2].legend(fontsize=8)

    for axis in axes:
        axis.grid(color="#d1d5db", linewidth=0.7, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)

    figure.suptitle(title, fontsize=15)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)

    for label, summary in summaries:
        print(
            f"{label}: duplicate_points={1.0 - summary['unique']:.6%} "
            f"zero_edges={summary['zero_edges']:.6%} "
            f"straight_points={summary['straight']:.6%} "
            f"center_points={summary['center']:.6%} "
            f"boundary_points={summary['boundary']:.6%} "
            f"median_nonzero_edge={summary['median_edge']:.10f}"
        )
    print(f"output={output}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare point placement in MA-CETSP populations"
    )
    parser.add_argument("instance", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("runs", nargs="+", type=parse_run)
    parser.add_argument("--title", default="MA-CETSP point placement")
    args = parser.parse_args()
    plot(read_instance(args.instance), args.runs, args.output, args.title)


if __name__ == "__main__":
    main()
