# Single-frame variant of generate.py, useful for quick checks.
#
# Reproducibility: same deterministic seeding scheme as generate.py
#   seed = 1_000_000 * L + 1000 * p_index + n
# The seed is set with np.random.seed() *inside* the numba nopython function,
# since numba keeps its own RNG state independent of NumPy's Python-level RNG.

import numpy as np
from numba import jit
import os

PC_INDEX = 500  # reserved p_index slot for the critical point p_c = 0.3116
PC_PROBABILITY = 0.3116


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


if __name__ == "__main__":
    """Generate a single percolation frame at p_c for L = 20."""

    L = 20
    coords = generate_sites(PC_PROBABILITY, L, frame_seed(L, PC_INDEX))

    os.makedirs(f"./{L}", exist_ok=True)
    lattice = f'Lattice="{L}.0 0.0 0.0 0.0 {L}.0 0.0 0.0 0.0 {L}.0"'
    with open(f"./{L}/percolation_sites_{PC_PROBABILITY:.4f}.xyz", "w") as f:
        f.write(f"{len(coords)}\n")
        f.write(f"{lattice}\n")
        for x, y, z in coords:
            f.write(f"1 {x} {y} {z}\n")
