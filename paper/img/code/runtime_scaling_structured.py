import numpy as np
import matplotlib.pyplot as plt

# --- Solver-only timing data (ms) ---
n = np.array([
    128, 256, 512, 1024, 2048, 4096,
    8192, 16384, 32768, 65536, 131072,
    262144, 524288, 1048576
], dtype=float)

time_ms = np.array([
    3, 14, 36, 54, 177, 407,
    592, 2045, 4179, 5683, 13007,
    41731, 65274, 161745
], dtype=float)

# --- Compute n log n and ratio ---
nlogn = n * np.log2(n)
ratio = time_ms / nlogn
C = np.mean(ratio)

# --- Plot runtime and normalization ---
plt.figure(figsize=(12,5))

# Raw runtime scaling
plt.subplot(1,2,1)
plt.plot(n, time_ms, 'o-', label='Measured runtime (ms)')
plt.plot(n, C * nlogn, '--', label='O(n log n) trend (scaled)')
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
plt.savefig("../runtime_scaling_structured_input.png", dpi=300)
plt.show()
