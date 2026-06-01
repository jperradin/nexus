# This script generates percolation sites in a simple cubic lattice.
#
# Reproducibility: every frame is generated from a deterministic, hardcoded seed
#   seed = 1_000_000 * L + 1000 * p_index + n
# where L is the lattice size, p_index the probability step (0..100 for the
# sweep, PC_INDEX for the critical point), and n the frame index. Because the
# generation runs inside a numba nopython function, the seed MUST be set with
# np.random.seed() *inside* that function -- numba keeps its own RNG state that
# is independent of NumPy's Python-level RNG.
#
# Each lattice size L also gets one serialized artifact ./{L}/percolation_sites_L{L}.npz
# holding every generated frame as an (M, 3) int32 array of occupied-site
# coordinates, so the exact dataset can be reloaded without regeneration.

import numpy as np
from numba import jit
import os
from tqdm import tqdm

PC_INDEX = 500  # reserved p_index slot for the critical point p_c = 0.3116
PC_PROBABILITY = 0.3116
N_FRAMES = 100


def frame_seed(L: int, p_index: int, n: int) -> int:
    """Deterministic, unique-per-frame seed. Stays below 2**32."""
    return 1_000_000 * L + 1000 * p_index + n


@jit(nopython=True, cache=True, fastmath=True)
def generate_sites(probability: float, size: int, seed: int) -> np.ndarray:
    np.random.seed(seed)
    coords = np.empty((size * size * size, 3), dtype=np.int32)
    count = 0
    for i in range(size):
        for j in range(size):
            for k in range(size):
                this_probability = np.random.random()
                if this_probability < probability:
                    coords[count, 0] = i % size
                    coords[count, 1] = j % size
                    coords[count, 2] = k % size
                    count += 1
    return coords[:count]


def write_frames_xyz(path: str, L: int, frames: list) -> None:
    """Write a multi-frame extended-XYZ trajectory (overwrite, not append)."""
    lattice = f'Lattice="{L}.0 0.0 0.0 0.0 {L}.0 0.0 0.0 0.0 {L}.0"'
    with open(path, "w") as f:
        for coords in frames:
            f.write(f"{len(coords)}\n")
            f.write(f"{lattice}\n")
            for x, y, z in coords:
                f.write(f"1 {x} {y} {z}\n")


if __name__ == "__main__":
    """
    Generate standard percolation lattice data sets for benchmarking validation:
        - At theoretical pc = 0.3116 for sizes from 20 to 50 with step 5
        - From 0.2 to 0.4 with step 0.002 for sizes from 20 to 50 with step 5
    """

    for L in tqdm(range(20, 55, 5)):
        os.makedirs(f"./{L}", exist_ok=True)
        archive = {}  # npz key -> (M, 3) int32 coordinates

        # Critical point: N_FRAMES independent frames at p = p_c
        pc_frames = []
        for n in range(N_FRAMES):
            coords = generate_sites(PC_PROBABILITY, L, frame_seed(L, PC_INDEX, n))
            pc_frames.append(coords)
            archive[f"pc_n{n:03d}"] = coords
        write_frames_xyz(
            f"./{L}/percolation_sites_{PC_PROBABILITY:.4f}.xyz", L, pc_frames
        )

        # Probability sweep for finite-size scaling
        for p_index in tqdm(range(0, 101), leave=False):
            probability = (p_index * 0.002) + 0.2
            sweep_frames = []
            for n in range(N_FRAMES):
                coords = generate_sites(probability, L, frame_seed(L, p_index, n))
                sweep_frames.append(coords)
                archive[f"p{p_index:03d}_n{n:03d}"] = coords
            write_frames_xyz(
                f"./{L}/percolation_sites_{probability:.4f}.xyz", L, sweep_frames
            )
