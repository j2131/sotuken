#!/usr/bin/env python3
# /// script
# dependencies = [
#   "numpy",
#   "matplotlib",
# ]
# ///
import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

def load_dat(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.read().strip().splitlines()
    if len(lines) < 2:
        return None, None
    time = float(lines[0].strip())
    nbody = int(lines[1].strip())
    data = []
    for line in lines[2:2 + nbody]:
        cols = line.strip().split()
        if len(cols) < 6:
            continue
        # id mass rad x y z ...
        data.append([float(cols[2]), float(cols[3]), float(cols[4]), float(cols[5])])
    if not data:
        return time, None
    arr = np.array(data)
    rad = arr[:, 0]
    pos = arr[:, 1:4]
    return time, (pos, rad)

def main():
    parser = argparse.ArgumentParser(description="Visualize DEM results (*.dat)")
    parser.add_argument("result_dir", nargs="?", default="result", help="result directory")
    parser.add_argument("--fps", type=int, default=30, help="animation fps")
    parser.add_argument("--index", type=int, default=None, help="show single frame index")
    parser.add_argument("--save", type=str, default=None, help="save animation to file (e.g. out.mp4 or out.gif)")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.result_dir, "*.dat")))
    if not files:
        raise SystemExit(f"No .dat files found in {args.result_dir}")

    if args.index is not None:
        idx = max(0, min(args.index, len(files) - 1))
        t, payload = load_dat(files[idx])
        if payload is None:
            raise SystemExit("No particle data")
        pos, rad = payload
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")
        # Scatter particles with their actual radius (scaling factor 1000 for visibility)
        ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], s=(rad * 1000) ** 2, alpha=0.7)
        ax.set_title(f"t={t:.4f}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        # Set bounds to see the whole cylinder area (R=1.0)
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.2, 1.2)
        ax.set_zlim(-1.2, 1.2)
        if args.save:
            plt.savefig(args.save)
            print(f"Saved image to {args.save}")
        else:
            plt.show()
        return

    frames = []
    times = []
    for path in files:
        t, payload = load_dat(path)
        if payload is None:
            continue
        pos, rad = payload
        frames.append((pos, rad))
        times.append(t)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    scat = ax.scatter([], [], [], s=[], alpha=0.7)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    # Set bounds to see the whole cylinder area
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(-1.2, 1.2)

    def init():
        scat._offsets3d = ([], [], [])
        scat.set_sizes([])
        return scat,

    def update(i):
        pos, rad = frames[i]
        scat._offsets3d = (pos[:, 0], pos[:, 1], pos[:, 2])
        scat.set_sizes((rad * 1000) ** 2)
        ax.set_title(f"t={times[i]:.4f}")
        return scat,

    ani = animation.FuncAnimation(
        fig, update, frames=len(frames), init_func=init, interval=1000 / args.fps, blit=False
    )

    if args.save:
        print(f"Saving animation to {args.save}...")
        if args.save.endswith(".gif"):
            ani.save(args.save, writer='pillow', fps=args.fps)
        else:
            ani.save(args.save, writer='ffmpeg', fps=args.fps)
        print("Done.")
    else:
        plt.show()

if __name__ == "__main__":
    main()
