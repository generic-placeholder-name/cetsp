import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

def draw_benchmark_shapes(input_file, output_file, image_dir):
    os.makedirs(image_dir, exist_ok=True)

    # Read circles from input file
    circles = []
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.strip():
                values = list(map(float, line.split()))
                if len(values) == 3:  # x, y, r
                    circles.append(values)

    # Read benchmark tour data
    tour_points = []
    total_length = 0.0

    with open(output_file, 'r') as outfile:
        lines = [line.strip() for line in outfile if line.strip()]
        if len(lines) >= 1:
            # Line 1: total length after ":"
            try:
                total_length = float(lines[0].split(":")[1].strip())
            except (IndexError, ValueError):
                print(f"Warning: Could not parse total length from {output_file}. Defaulting to 0.")
                total_length = 0.0

            # Lines 3+ are node, x, y
            for line in lines[2:]:
                parts = line.split()
                if len(parts) >= 3:
                    x, y = float(parts[1]), float(parts[2])
                    tour_points.append([x, y])

    # Plotting
    fig, ax = plt.subplots()
    min_x, min_y = float('inf'), float('inf')
    max_x, max_y = float('-inf'), float('-inf')

    # Draw circles
    for x, y, r in circles:
        circle = patches.Circle((x, y), r, edgecolor='blue', facecolor='none', linewidth=0.5)
        ax.add_patch(circle)
        min_x, min_y = min(min_x, x - r), min(min_y, y - r)
        max_x, max_y = max(max_x, x + r), max(max_y, y + r)

    # Draw tour
    if tour_points:
        for i in range(len(tour_points)):
            x1, y1 = tour_points[i]
            x2, y2 = tour_points[(i + 1) % len(tour_points)]
            ax.plot([x1, x2], [y1, y2], 'r-', linewidth=1)

    ax.set_aspect('equal', 'box')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(f'Tour and Circles\nTotal Tour Length: {total_length:.2f}')
    ax.grid(True)
    ax.set_xlim(min_x - 1, max_x + 1)
    ax.set_ylim(min_y - 1, max_y + 1)

    output_image = os.path.join(image_dir, f'{os.path.basename(input_file)}_plot.png')
    plt.savefig(output_image)
    plt.close(fig)
    print(f"Plot saved to {output_image}")


# Main processing
input_dir = './cmu_processed'    # Same as before
output_dir = './benchmark'       # Benchmark output files
image_dir = './benchmark_images'

for filename in os.listdir(input_dir):
    if filename.endswith('.txt'):
        input_file = os.path.join(input_dir, filename)
        output_file = os.path.join(output_dir, filename)
        if os.path.exists(output_file):
            draw_benchmark_shapes(input_file, output_file, image_dir)
