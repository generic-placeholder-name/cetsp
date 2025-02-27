import os
import re

# Define the directory containing the text files
input_dir = './cmu'
output_dir = './cmu_processed'

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Iterate over all files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith('.cetsp'):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename.replace(".cetsp", ".txt"))

        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                # Skip lines starting with '//'
                if line.strip().startswith('//'):
                    # Check if the line starts with //Depot and contains numbers
                    if line.strip().startswith('//Depot'):
                        # Find numbers after the "//Depot" part
                        match = re.search(r'[-+]?\d*\.\d+|\d+', line)
                        if match:
                            # Extract the first three numbers from the line
                            numbers = re.findall(r'[-+]?\d*\.\d+|\d+', line)
                            if len(numbers) >= 3:
                                x, y, z = numbers[0], numbers[1], 0  # Set z to 0
                                outfile.write(f'{x} {y} {z}\n')
                    continue  # Skip the line if it's a comment or processed line

                # Skip empty lines
                if line.strip() == '':
                    continue

                # Check if the line contains numbers
                parts = line.split()
                try:
                    # Extract the numerical values from the line
                    numbers = [float(part) for part in parts]
                    # Select a, b, and d (1st, 2nd, and 4th numbers)
                    filtered_numbers = [str(numbers[0]), str(numbers[1]), str(numbers[3])]
                    outfile.write(' '.join(filtered_numbers) + '\n')
                except ValueError:
                    # If line doesn't contain valid numbers, skip it
                    pass
