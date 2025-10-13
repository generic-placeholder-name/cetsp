import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

# Function to read the files and draw the plot
def draw_shapes_from_files(input_file, output_file, image_dir):
    # Create image directory if it doesn't exist
    os.makedirs(image_dir, exist_ok=True)

    # Read the circle data from the input file
    circles = []
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.strip():  # Skip empty lines
                values = list(map(float, line.split()))
                if len(values) == 3:  # Ensure the line has 3 numbers (x, y, r)
                    circles.append(values)

    # Read the tour data and total distance from the output file
    tour_points = []
    total_length = 0.0

    with open(output_file, 'r') as outfile:
        lines = [line.strip() for line in outfile if line.strip()]  # Skip empty lines
        if lines:
            try:
                total_length = float(lines[0])  # First line = total tour length
            except ValueError:
                print(f"Warning: Could not parse total length from {output_file}. Defaulting to 0.")
                total_length = 0.0

            # Remaining lines are tour points
            for line in lines[1:]:
                values = list(map(float, line.split()))
                if len(values) == 2:
                    tour_points.append(values)

    # Prepare to plot
    fig, ax = plt.subplots()

    # Variables to track min/max for the plot limits
    min_x, min_y = float('inf'), float('inf')
    max_x, max_y = float('-inf'), float('-inf')

    # Draw circles
    for (x, y, r) in circles:
        circle = patches.Circle((x, y), r, edgecolor='blue', facecolor='none', linewidth=0.5)
        ax.add_patch(circle)

        # Update min/max for axis limits based on the circle centers and radii
        min_x = min(min_x, x - r)
        min_y = min(min_y, y - r)
        max_x = max(max_x, x + r)
        max_y = max(max_y, y + r)

    # Draw the tour
    if tour_points:
        for i in range(len(tour_points)):
            x1, y1 = tour_points[i]
            x2, y2 = tour_points[(i + 1) % len(tour_points)]  # Wrap around to the first point
            ax.plot([x1, x2], [y1, y2], 'r-', linewidth=1)

    # Set the aspect ratio to be equal
    ax.set_aspect('equal', 'box')

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(f'Tour and Circles\nTotal Tour Length: {total_length:.2f}')

    # Set a grid
    ax.grid(True)

    # Set limits for the axes based on the circle centers and their radii
    ax.set_xlim(min_x - 1, max_x + 1)
    ax.set_ylim(min_y - 1, max_y + 1)

    # Save the plot to a file
    output_image = os.path.join(image_dir, f'{os.path.basename(input_file)}_plot.png')
    plt.savefig(output_image)

    # Close the plot
    plt.close(fig)
    print(f"Plot saved to {output_image}")


# Process all files in the input and output directories
input_dir = './cmu_processed'  # Input directory with circle files
output_dir = './cmu_output'    # Output directory with tour files
image_dir = './output_images'  # Directory for the generated images

# Iterate over all files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith('.txt'):  # Only process .txt files
        input_file = os.path.join(input_dir, filename)
        output_file = os.path.join(output_dir, filename)  # Corresponding file in output_dir
        if os.path.exists(output_file):  # Ensure the corresponding tour file exists
            draw_shapes_from_files(input_file, output_file, image_dir)
