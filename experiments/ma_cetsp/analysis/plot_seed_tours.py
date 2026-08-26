#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
from statistics import mean

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize


def read_instance(path):
    tokens = path.read_text(encoding="utf-8").split()
    if not tokens:
        raise ValueError(f"empty instance file: {path}")

    node_count = int(tokens[0])
    if len(tokens) != 1 + 3 * node_count:
        raise ValueError(f"invalid instance file: {path}")

    circles = []
    offset = 1
    for _ in range(node_count):
        circles.append(tuple(float(value) for value in tokens[offset : offset + 3]))
        offset += 3
    return circles


def read_seeds(path):
    tokens = path.read_text(encoding="utf-8").split()
    if len(tokens) < 4 or tokens[:2] != ["MA_CETSP_SEEDS", "1"]:
        raise ValueError(f"invalid seed file header: {path}")

    node_count = int(tokens[2])
    tour_count = int(tokens[3])
    tours = []
    offset = 4
    for _ in range(tour_count):
        if offset >= len(tokens) or tokens[offset] != "TOUR":
            raise ValueError(f"missing TOUR record in {path}")
        offset += 1

        tour = []
        for _ in range(node_count):
            if offset + 2 >= len(tokens):
                raise ValueError(f"truncated TOUR record in {path}")
            tour.append(
                (int(tokens[offset]), float(tokens[offset + 1]), float(tokens[offset + 2]))
            )
            offset += 3
        tours.append(tour)

    if offset != len(tokens):
        raise ValueError(f"unexpected data after final TOUR record in {path}")
    return node_count, tours


def closed_points(tour):
    points = [(x, y) for _, x, y in tour]
    return points + points[:1]


def tour_length(tour):
    points = closed_points(tour)
    return sum(
        math.hypot(x1 - x0, y1 - y0)
        for (x0, y0), (x1, y1) in zip(points, points[1:])
    )


def tour_edges(tour):
    ids = [node_id for node_id, _, _ in tour]
    return {
        tuple(sorted((first, second)))
        for first, second in zip(ids, ids[1:] + ids[:1])
    }


def percentile(values, fraction):
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def diversity_summary(circles, tours):
    edge_sets = [tour_edges(tour) for tour in tours]
    pairwise_overlap = []
    for first in range(len(edge_sets)):
        for second in range(first + 1, len(edge_sets)):
            pairwise_overlap.append(
                len(edge_sets[first] & edge_sets[second]) / len(edge_sets[first])
            )

    all_edges = set().union(*edge_sets)
    edge_lengths = []
    for first, second in all_edges:
        x0, y0, _ = circles[first]
        x1, y1, _ = circles[second]
        edge_lengths.append(math.hypot(x1 - x0, y1 - y0))
    long_threshold = percentile(edge_lengths, 0.75)

    long_edge_sets = []
    for edges in edge_sets:
        long_edges = set()
        for edge in edges:
            first, second = edge
            x0, y0, _ = circles[first]
            x1, y1, _ = circles[second]
            if math.hypot(x1 - x0, y1 - y0) >= long_threshold:
                long_edges.add(edge)
        long_edge_sets.append(long_edges)

    long_jaccard = []
    for first in range(len(long_edge_sets)):
        for second in range(first + 1, len(long_edge_sets)):
            union = long_edge_sets[first] | long_edge_sets[second]
            if union:
                long_jaccard.append(
                    len(long_edge_sets[first] & long_edge_sets[second]) / len(union)
                )

    return {
        "distinct": len({tuple(node_id for node_id, _, _ in tour) for tour in tours}),
        "edge_overlap_mean": mean(pairwise_overlap),
        "edge_overlap_min": min(pairwise_overlap),
        "edge_overlap_max": max(pairwise_overlap),
        "long_threshold": long_threshold,
        "long_jaccard_mean": mean(long_jaccard),
    }


def set_map_style(ax, bounds):
    min_x, max_x, min_y, max_y = bounds
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#d1d5db")


def draw_background(ax, circles):
    ax.scatter(
        [circle[0] for circle in circles],
        [circle[1] for circle in circles],
        s=2.0,
        color="#9ca3af",
        alpha=0.38,
        linewidths=0,
        zorder=1,
    )
    ax.scatter(
        [circles[0][0]],
        [circles[0][1]],
        marker="*",
        s=34,
        color="#dc2626",
        edgecolors="white",
        linewidths=0.35,
        zorder=4,
    )


def plot(circles, tours, output, title):
    columns = min(5, len(tours))
    tour_rows = math.ceil(len(tours) / columns)
    figure = plt.figure(
        figsize=(3.45 * columns, 3.25 * tour_rows + 5.0),
        constrained_layout=True,
    )
    grid = figure.add_gridspec(
        tour_rows + 1,
        columns,
        height_ratios=[1] * tour_rows + [1.4],
    )

    min_x = min(x - radius for x, _, radius in circles)
    max_x = max(x + radius for x, _, radius in circles)
    min_y = min(y - radius for _, y, radius in circles)
    max_y = max(y + radius for _, y, radius in circles)
    padding = 0.025 * max(max_x - min_x, max_y - min_y)
    bounds = (min_x - padding, max_x + padding, min_y - padding, max_y + padding)

    summary = diversity_summary(circles, tours)
    lengths = [tour_length(tour) for tour in tours]
    colors = plt.colormaps["viridis"]

    for index, tour in enumerate(tours):
        row, column = divmod(index, columns)
        ax = figure.add_subplot(grid[row, column])
        draw_background(ax, circles)
        points = closed_points(tour)
        ax.plot(
            [point[0] for point in points],
            [point[1] for point in points],
            color=colors(index / max(1, len(tours) - 1)),
            linewidth=0.72,
            alpha=0.92,
            zorder=2,
        )
        ax.set_title(f"Tour {index + 1}  |  {lengths[index]:.2f}", fontsize=9)
        set_map_style(ax, bounds)

    for index in range(len(tours), tour_rows * columns):
        row, column = divmod(index, columns)
        figure.add_subplot(grid[row, column]).axis("off")

    overlay_columns = max(1, columns // 2)
    overlay_ax = figure.add_subplot(grid[tour_rows, :overlay_columns])
    draw_background(overlay_ax, circles)
    for tour in tours:
        points = closed_points(tour)
        overlay_ax.plot(
            [point[0] for point in points],
            [point[1] for point in points],
            color="#2563eb",
            linewidth=0.65,
            alpha=0.075,
            zorder=2,
        )
    overlay_ax.set_title("All optimized paths (darker = repeated geometry)", fontsize=10)
    set_map_style(overlay_ax, bounds)

    consensus_ax = figure.add_subplot(grid[tour_rows, overlay_columns:])
    draw_background(consensus_ax, circles)
    edge_counts = {}
    for tour in tours:
        for edge in tour_edges(tour):
            edge_counts[edge] = edge_counts.get(edge, 0) + 1

    ordered_edges = sorted(edge_counts, key=edge_counts.get)
    segments = []
    counts = []
    widths = []
    for first, second in ordered_edges:
        segments.append(
            [(circles[first][0], circles[first][1]), (circles[second][0], circles[second][1])]
        )
        count = edge_counts[(first, second)]
        counts.append(count)
        widths.append(0.2 + 2.2 * count / len(tours))

    collection = LineCollection(
        segments,
        array=counts,
        cmap="plasma",
        norm=Normalize(1, len(tours)),
        linewidths=widths,
        alpha=0.72,
        zorder=2,
    )
    consensus_ax.add_collection(collection)
    consensus_ax.set_title("Ordering consensus on circle centers", fontsize=10)
    set_map_style(consensus_ax, bounds)
    colorbar = figure.colorbar(collection, ax=consensus_ax, fraction=0.035, pad=0.02)
    colorbar.set_label("Tours sharing edge", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)

    figure.suptitle(
        f"{title}\n"
        f"{summary['distinct']}/{len(tours)} distinct orders; mean pairwise shared edges "
        f"{summary['edge_overlap_mean']:.1%} "
        f"(range {summary['edge_overlap_min']:.1%}–{summary['edge_overlap_max']:.1%}); "
        f"long-edge Jaccard {summary['long_jaccard_mean']:.1%} "
        f"(center distance ≥ {summary['long_threshold']:.1f})",
        fontsize=13,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=170, bbox_inches="tight")
    plt.close(figure)
    return summary, lengths


def main():
    parser = argparse.ArgumentParser(description="Plot every tour in an MA-CETSP seed population")
    parser.add_argument("instance", type=Path)
    parser.add_argument("seeds", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--title", default="MA-CETSP seed population")
    args = parser.parse_args()

    circles = read_instance(args.instance)
    node_count, tours = read_seeds(args.seeds)
    if node_count != len(circles):
        raise ValueError(
            f"seed node count {node_count} does not match instance count {len(circles)}"
        )

    summary, lengths = plot(circles, tours, args.output, args.title)
    print(f"output={args.output}")
    print(f"tours={len(tours)}")
    print(f"distinct_orders={summary['distinct']}")
    print(f"best_length={min(lengths):.10f}")
    print(f"worst_length={max(lengths):.10f}")
    print(f"mean_pairwise_shared_edges={summary['edge_overlap_mean']:.10f}")
    print(f"long_edge_threshold={summary['long_threshold']:.10f}")
    print(f"mean_pairwise_long_edge_jaccard={summary['long_jaccard_mean']:.10f}")


if __name__ == "__main__":
    main()
