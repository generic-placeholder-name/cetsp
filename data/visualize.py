import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

# Function to read the file and draw the plot
def draw_shapes_from_file(input_file, output_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read the data from the file
    points = []
    total_length = None  # Initialize total_length variable

    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        
        # Skip the first line (header or irrelevant info)
        for line in lines[1:-1]:  # Process all lines except the first and the last one
            if line.strip():  # Skip empty lines
                values = list(map(float, line.split()))
                if len(values) == 5:  # Ensure the line has 5 numbers
                    points.append(values)

        # Extract the total length from the last line
        last_line = lines[-1].strip()
        if last_line.startswith("Total tour distance:"):
            try:
                total_length = float(last_line.split(":")[1].strip())  # Extract number after the colon
            except ValueError:
                print("The last line does not contain a valid number for total length.")
    
    # Prepare to plot
    fig, ax = plt.subplots()

    # List to store the second points for the lines
    line_points = []

    # Variables to track min/max for the plot limits
    min_x, min_y = float('inf'), float('inf')
    max_x, max_y = float('-inf'), float('-inf')

    # Draw circles, lines, and points
    for (x1, y1, d, x2, y2) in points:
        # Draw the circle centered at (x1, y1) with radius d
        circle = patches.Circle((x1, y1), d, edgecolor='blue', facecolor='none', linewidth=0.5)
        ax.add_patch(circle)

        # Plot the black point at (x1, y1)
        # ax.plot(x1, y1, 'ko', markersize=6)

        # Plot the green point at (x2, y2)
        # ax.plot(x2, y2, 'go', markersize=6)

        # Append (x2, y2) to the line_points list for connecting lines
        line_points.append((x2, y2))

        # Update min/max for axis limits based on the circle centers and radii
        min_x = min(min_x, x1 - d)
        min_y = min(min_y, y1 - d)
        max_x = max(max_x, x1 + d)
        max_y = max(max_y, y1 + d)

    # Connect the points in line_points list with red lines (connect second points)
    if line_points:
        for i in range(len(line_points) - 1):
            x1, y1 = line_points[i]
            x2, y2 = line_points[i + 1]
            ax.plot([x1, x2], [y1, y2], 'r-', linewidth=1)

        # Connect the last point back to the first one
        x1, y1 = line_points[-1]
        x2, y2 = line_points[0]
        ax.plot([x1, x2], [y1, y2], 'r-', linewidth=1)

    # Set the aspect ratio to be equal
    ax.set_aspect('equal', 'box')

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(f'Circles, Lines, and Points\nTotal Tour Length: {total_length}' if total_length else 'Circles, Lines, and Points')

    # Set a grid
    ax.grid(True)

    # Set limits for the axes based on the circle centers and their radii
    ax.set_xlim(min_x - 1, max_x + 1)
    ax.set_ylim(min_y - 1, max_y + 1)

    # Save the plot to a file
    output_file = os.path.join(output_dir, f'{os.path.basename(input_file)}_plot.png')
    plt.savefig(output_file)

    # Close the plot
    plt.close(fig)
    print(f"Plot saved to {output_file}")

# Process all files in the input directory
input_dir = './cmu_output'  # Input directory with the files
output_dir = './output_images'   # Output directory for the images

# Iterate over all files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith('.txt'):  # Only process .txt files
        input_file = os.path.join(input_dir, filename)
        draw_shapes_from_file(input_file, output_dir)
