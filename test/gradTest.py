import numpy as np
import matplotlib.pyplot as plt
import math

# --- Problem Generation ---
def generate_problem_instance(num_circles=10, world_size=100):
    core_radius = world_size * 0.1
    core_center = np.array([0.0, 0.0])
    circles = []
    for _ in range(num_circles):
        offset = (np.random.rand(2) - 0.5) * world_size * 0.6
        center = core_center + offset
        dist = np.linalg.norm(center - core_center)
        min_r = dist + core_radius
        r = min_r + np.random.rand() * world_size * 0.2
        circles.append((center, r))
    p0 = core_center.copy()
    def seg_dist(p, a, b):
        ab = b - a
        t = np.dot(p - a, ab) / np.dot(ab, ab)
        t = np.clip(t, 0.0, 1.0)
        proj = a + t * ab
        return np.linalg.norm(p - proj)
    while True:
        θ1, θ2 = np.random.rand() * 2*np.pi, np.random.rand() * 2*np.pi
        r1 = world_size * (0.6 + np.random.rand() * 0.4)
        r2 = world_size * (0.6 + np.random.rand() * 0.4)
        a = np.array([np.cos(θ1), np.sin(θ1)]) * r1
        b = np.array([np.cos(θ2), np.sin(θ2)]) * r2
        if seg_dist(core_center, a, b) > core_radius:
            break
    return circles, p0, a, b

# --- Gradient Computation ---
def gradient_and_max_step(p, a, b, circles):
    va, vb = a - p, b - p
    da, db = np.linalg.norm(va), np.linalg.norm(vb)
    if da < 1e-9 or db < 1e-9:
        return None, 0.0
    g = va/da + vb/db
    gn = np.linalg.norm(g)
    if gn < 1e-9:
        return None, 0.0
    d = g / gn
    t_max = np.inf
    for c, r in circles:
        pc = p - c
        B = 2 * np.dot(d, pc)
        C = np.dot(pc, pc) - r*r
        disc = B*B - 4*C
        if disc >= 0:
            sd = math.sqrt(disc)
            t1, t2 = (-B + sd)/2, (-B - sd)/2
            t_pos = max(t1, t2)
            t_max = min(t_max, t_pos)
    return d, t_max

# --- Descent Methods ---
def descent_max(circles, p0, a, b, steps):
    p, path = p0.copy(), [p0.copy()]
    for _ in range(steps):
        d, t = gradient_and_max_step(p, a, b, circles)
        if d is None or t <= 0:
            break
        p = p + t * d
        path.append(p.copy())
    return np.array(path)

def descent_standard(circles, p0, a, b, steps, lr):
    p, path = p0.copy(), [p0.copy()]
    for _ in range(steps):
        d, t = gradient_and_max_step(p, a, b, circles)
        if d is None or t <= 0:
            break
        step = min(lr, t)
        p = p + step * d
        path.append(p.copy())
    return np.array(path)

def tune_lr(circles, p0, a, b, steps, candidates):
    best_lr, best_dist = None, np.inf
    for lr in candidates:
        path = descent_standard(circles, p0, a, b, steps, lr)
        final = path[-1]
        dist = np.linalg.norm(a-final) + np.linalg.norm(b-final)
        if dist < best_dist:
            best_dist, best_lr = dist, lr
    print(f"Tuned lr: {best_lr:.3f} (final distance {best_dist:.4f})")
    return best_lr

def descent_half_then_max(circles, p0, a, b, k):
    p, path = p0.copy(), [p0.copy()]
    for _ in range(k):
        d, t = gradient_and_max_step(p, a, b, circles)
        if d is None or t <= 0:
            break
        p = p + 0.5 * t * d
        path.append(p.copy())
    d, t = gradient_and_max_step(p, a, b, circles)
    if d is not None and t > 0:
        p = p + t * d
        path.append(p.copy())
    return np.array(path)

# --- Plotting and Saving ---
def plot_and_save(circles, a, b, paths, labels, prefix="plot"):
    # Paths figure
    fig1, ax1 = plt.subplots(figsize=(8,8))
    for c, r in circles:
        ax1.add_patch(plt.Circle(c, r, color='skyblue', alpha=0.3))
    ax1.plot(a[0], a[1], 'go', label='a')
    ax1.plot(b[0], b[1], 'ro', label='b')
    for path, lab in zip(paths, labels):
        ax1.plot(path[:,0], path[:,1], '-o', label=lab)
    ax1.set_aspect('equal'); ax1.legend(); ax1.grid(True)
    ax1.set_title('Paths of Descent Methods')
    file1 = f"{prefix}_paths.png"
    fig1.savefig(file1)
    print(f"Saved paths plot to {file1}")

    # Distance vs Iteration figure
    fig2, ax2 = plt.subplots(figsize=(8,6))
    dist_lists = []
    for path in paths:
        dist_lists.append([np.linalg.norm(a-p)+np.linalg.norm(b-p) for p in path])
    # build legend labels with final distances
    legend_labels = [f"{lab} (final {dists[-1]:.2f})" for lab, dists in zip(labels, dist_lists)]
    for dists, lab in zip(dist_lists, legend_labels):
        ax2.plot(dists, '-o', label=lab)
    ax2.set_xlabel('Iteration'); ax2.set_ylabel('Distance'); ax2.legend(); ax2.grid(True)
    ax2.set_title('Distance vs Iteration')
    file2 = f"{prefix}_distance.png"
    fig2.savefig(file2)
    print(f"Saved distance plot to {file2}")

# --- Main ---
def main():
    NUM_CIRCLES, STEPS, WORLD = 10, 5, 100
    circles, p0, a, b = generate_problem_instance(NUM_CIRCLES, WORLD)

    # Compute paths
    p_max = descent_max(circles, p0, a, b, STEPS)
    lr_candidates = np.logspace(-2, 1, 8)
    best_lr = tune_lr(circles, p0, a, b, STEPS, lr_candidates)
    p_std = descent_standard(circles, p0, a, b, STEPS, best_lr)
    K = STEPS - 1
    p_hyb = descent_half_then_max(circles, p0, a, b, K)

    # Plot & save
    plot_and_save(circles, a, b,
                  [p_max, p_std, p_hyb],
                  ['Max-Step', f'Std lr={best_lr:.2f}', f'Half({K})+Max'],
                  prefix="steiner_descent")

if __name__ == '__main__':
    main()
