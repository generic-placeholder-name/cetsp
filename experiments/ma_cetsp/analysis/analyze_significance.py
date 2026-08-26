#!/usr/bin/env python3

import argparse
import csv
import math
import statistics
from pathlib import Path

from scipy.stats import t as student_t


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare the two 20-run Mennell protocols with Welch tests."
    )
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
        "--output",
        type=Path,
        default=Path("results/mennell_hybrid/tables/mennell_protocol_significance.csv"),
    )
    return parser.parse_args()


def welch_test(mean_a, sd_a, n_a, mean_b, sd_b, n_b):
    variance_a = sd_a * sd_a / n_a
    variance_b = sd_b * sd_b / n_b
    standard_error_squared = variance_a + variance_b
    if standard_error_squared == 0.0:
        return (0.0, math.inf, 1.0) if mean_a == mean_b else (
            math.copysign(math.inf, mean_a - mean_b),
            math.inf,
            0.0,
        )

    statistic = (mean_a - mean_b) / math.sqrt(standard_error_squared)
    denominator = 0.0
    if variance_a:
        denominator += variance_a * variance_a / (n_a - 1)
    if variance_b:
        denominator += variance_b * variance_b / (n_b - 1)
    degrees_of_freedom = standard_error_squared * standard_error_squared / denominator
    p_value = 2.0 * student_t.sf(abs(statistic), degrees_of_freedom)
    return statistic, degrees_of_freedom, p_value


def holm_adjust(rows):
    ordered = sorted(rows, key=lambda row: row["p_value"])
    previous = 0.0
    count = len(ordered)
    for index, row in enumerate(ordered):
        adjusted = min(1.0, (count - index) * row["p_value"])
        row["holm_adjusted_p_value"] = max(previous, adjusted)
        previous = row["holm_adjusted_p_value"]


def main():
    args = parse_args()
    runs_by_instance = {}
    with args.runs.open(newline="", encoding="utf-8-sig") as input_file:
        for row in csv.DictReader(input_file):
            runs_by_instance.setdefault(row["instance"], []).append(float(row["final"]))

    rows = []
    with args.lei_hao.open(newline="", encoding="utf-8-sig") as input_file:
        for published in csv.DictReader(input_file):
            instance = published["instance"]
            if instance not in runs_by_instance:
                raise ValueError(f"No hybrid runs for {instance}")
            ours = runs_by_instance.pop(instance)
            if len(ours) != 20:
                raise ValueError(f"Expected 20 hybrid runs for {instance}, got {len(ours)}")
            ours_mean = statistics.mean(ours)
            ours_sd = statistics.stdev(ours)
            lei_hao_mean = float(published["lei_hao_mean"])
            lei_hao_sd = float(published["lei_hao_sd"])
            lei_hao_runs = int(published["lei_hao_runs"])
            # Tables 5--7 report both statistics to two decimal places. Use
            # the point in each rounding interval that is least favorable to
            # rejecting equal means.
            rounding_half_width = 0.005
            reference_mean = min(
                max(ours_mean, lei_hao_mean - rounding_half_width),
                lei_hao_mean + rounding_half_width,
            )
            reference_sd = lei_hao_sd + rounding_half_width
            statistic, degrees_of_freedom, p_value = welch_test(
                ours_mean,
                ours_sd,
                len(ours),
                reference_mean,
                reference_sd,
                lei_hao_runs,
            )
            rows.append(
                {
                    "instance": instance,
                    "our_mean": ours_mean,
                    "our_sd": ours_sd,
                    "our_runs": len(ours),
                    "lei_hao_mean": lei_hao_mean,
                    "lei_hao_sd": lei_hao_sd,
                    "lei_hao_runs": lei_hao_runs,
                    "mean_difference": ours_mean - lei_hao_mean,
                    "rounding_aware_mean_difference": ours_mean - reference_mean,
                    "reference_mean_for_test": reference_mean,
                    "reference_sd_for_test": reference_sd,
                    "welch_t": statistic,
                    "degrees_of_freedom": degrees_of_freedom,
                    "p_value": p_value,
                    "holm_adjusted_p_value": 1.0,
                    "classification": "inconclusive",
                }
            )

    if runs_by_instance:
        raise ValueError(f"Missing published statistics for: {sorted(runs_by_instance)}")

    holm_adjust(rows)
    for row in rows:
        if row["holm_adjusted_p_value"] < 0.05:
            row["classification"] = (
                "significant_win" if row["mean_difference"] < 0.0 else "significant_loss"
            )

    fieldnames = list(rows[0])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["instance"]))

    counts = {classification: 0 for classification in (
        "significant_win", "significant_loss", "inconclusive"
    )}
    for row in rows:
        counts[row["classification"]] += 1
    print(", ".join(f"{key}={value}" for key, value in counts.items()))


if __name__ == "__main__":
    main()
