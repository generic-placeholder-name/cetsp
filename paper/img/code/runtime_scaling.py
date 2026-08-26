import numpy as np
import matplotlib.pyplot as plt

# --- Solver-only timing data (ms) ---
n = np.array([
    128, 256, 512, 1024, 2048, 4096,
    8192, 16384, 32768, 65536, 131072,
    262144, 524288, 1048576
], dtype=float)

time_ms = np.array([
    46, 116, 248, 550, 1264, 2505,
    5161, 9578, 16677, 30758, 49155,
    74206, 121820, 198618
], dtype=float)

# --- Compute n log n and ratio ---
nlogn = n * np.log2(n)
ratio = time_ms / nlogn

# --- Plot runtime and normalization ---
plt.figure(figsize=(12,5))

# Raw runtime scaling
plt.subplot(1,2,1)
plt.plot(n, time_ms, 'o-', label='Measured runtime (ms)')
plt.plot(n, np.mean(ratio) * nlogn, '--', label='O(n log n) trend (scaled)')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('n')
plt.ylabel('Runtime (ms)')
plt.title('Runtime scaling')
plt.legend()
plt.grid(which='major', linestyle='-', linewidth=0.5, alpha=0.5)   # major grid
plt.grid(which='minor', linestyle=':', linewidth=0.2, alpha=0.2)    # minor grid

# Normalized runtime
plt.subplot(1,2,2)
plt.plot(n, ratio, 'o-')
plt.xscale('log')
plt.xlabel('n')
plt.ylabel('time / (n log n)')
plt.title('Normalized runtime (divided by n log n)')
plt.grid(which='major', linestyle='-', linewidth=0.5, alpha=0.5)
plt.grid(which='minor', linestyle=':', linewidth=0.2, alpha=0.2)

plt.tight_layout()
plt.savefig("../runtime_scaling_random_input.png", dpi=300)
plt.show()
