import os
import numpy as np

input_dir = 'cmu'
benchmark_dir = 'benchmark'
output_dir = 'cmu_processed'
os.makedirs(output_dir, exist_ok=True)
files = ['d493.cetsp', 'dsj1000.cetsp', 'kroD100.cetsp', 'lin318.cetsp', 'pcb442.cetsp', 'rat195.cetsp', 'rd400.cetsp']
overlap_ratios = [0.02, 0.1, 0.3]

for fname in files:
    base_path = os.path.join(input_dir, fname)
    base_points = []
    depot = None

    # Read base file
    with open(base_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('//Depot:'):
                depot = tuple(map(float, line.split(':')[1].split(',')[:2]))
            elif not line.startswith('//'):
                parts = line.split(' ')
                if len(parts) >= 5:
                    x, y, z, r, id_ = parts
                    base_points.append((float(x), float(y)))

    if depot is not None:
        base_points.append(depot)
    base_points = np.array(base_points)
    N = len(base_points)
    if N == 0:
        print(f"Skipping {fname}: no points or depot not found.")
        continue

    for oratio in overlap_ratios:
        # Load tour points from benchmark
        tour_name = os.path.splitext(fname)[0] + f"_or{int(oratio*100)}.txt"
        tour_path = os.path.join(benchmark_dir, tour_name)
        if not os.path.exists(tour_path):
            print(f"Skipping {tour_name}: not found in benchmark_dir.")
            continue
        tour_points = []
        with open(tour_path, 'r') as tf:
            lines = tf.readlines()
            for line in lines[2:]:  # skip first two lines
                parts = line.strip().split()
                if len(parts) >= 3:
                    x, y = map(float, parts[1:3])
                    tour_points.append((x, y))
        tour_points = np.array(tour_points)
        if len(tour_points) == 0:
            print(f"Skipping {tour_name}: no tour points found.")
            continue

        # For each base point, find closest tour point
        dists = []
        for bp in base_points[:-1]:  # exclude depot
            dist = np.min(np.linalg.norm(tour_points - bp, axis=1))
            dists.append(dist)
        r_final = max(dists)

        # Write output
        outname = os.path.splitext(fname)[0] + f"_or{int(oratio*100)}.txt"
        outpath = os.path.join(output_dir, outname)
        with open(outpath, 'w') as out:
            for (x, y) in base_points[:-1]:  # all except depot
                out.write(f"{x} {y} {r_final:.6f}\n")
            # Write depot last
            x, y = base_points[-1]
            out.write(f"{x} {y} 0\n")
