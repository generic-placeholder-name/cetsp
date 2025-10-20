import os
import csv

BENCHMARK_DIR = "./benchmark"
OUTPUT_DIR = "./cmu_output"
OUTPUT_CSV = "comparison.csv"

def read_benchmark_value(filepath):
    with open(filepath, 'r') as f:
        line = f.readline().strip()
        # Expected format: "value : x"
        try:
            return float(line.split(":")[1].strip())
        except (IndexError, ValueError):
            print(f"Warning: Could not parse benchmark value from {filepath}")
            return None

def read_output_value(filepath):
    with open(filepath, 'r') as f:
        line = f.readline().strip()
        try:
            return float(line)
        except ValueError:
            print(f"Warning: Could not parse output value from {filepath}")
            return None

def main():
    benchmark_files = sorted(os.listdir(BENCHMARK_DIR))
    results = []

    for filename in benchmark_files:
        bench_path = os.path.join(BENCHMARK_DIR, filename)
        output_path = os.path.join(OUTPUT_DIR, filename)

        if not os.path.isfile(output_path):
            print(f"Warning: No corresponding output file for {filename}")
            continue

        bench_value = read_benchmark_value(bench_path)
        output_value = read_output_value(output_path)

        if bench_value is None or output_value is None:
            continue

        percent_diff = ((output_value - bench_value) / bench_value) * 100 if bench_value != 0 else float('inf')

        results.append({
            "filename": filename,
            "benchmark_value": bench_value,
            "output_value": output_value,
            "percent_diff": percent_diff
        })

    with open(OUTPUT_CSV, 'w', newline='') as csvfile:
        fieldnames = ["filename", "benchmark_value", "output_value", "percent_diff"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"Comparison complete. Results saved to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
