import os
import random
import math

def generate_cetsp_tests(output_dir="timing_input"):
    os.makedirs(output_dir, exist_ok=True)

    for exp in range(7, 21):  # 2^7 = 128, 2^20 = 1,048,576
        n = 2 ** exp
        filename = os.path.join(output_dir, f"test_n{n}.txt")

        with open(filename, "w") as f:
            for _ in range(n):
                # Random center in [-n, n]
                cx = random.uniform(-n, n)
                cy = random.uniform(-n, n)

                # Radius randomized based on n
                radius = n * random.uniform(0.01, 0.02)

                f.write(f"{cx:.6f} {cy:.6f} {radius:.6f}\n")

        print(f"Generated {filename} with {n} circles.")

if __name__ == "__main__":
    generate_cetsp_tests()
