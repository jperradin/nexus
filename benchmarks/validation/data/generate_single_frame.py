# This script generate percolation sites in a 100x100x100 orthohombic lattice

import numpy as np
from numba import jit
from typing import List
import os
from tqdm import tqdm


@jit(nopython=True, cache=True, fastmath=True)
def generate_sites(probability: float, size: int) -> List[str]:
    sites = []
    for i in range(size):
        for j in range(size):
            for k in range(size):
                this_probability = np.random.random()
                # generate the lattice position of the sites
                x = i % size
                y = j % size
                z = k % size

                if this_probability < probability:
                    sites.append(f"1 {x} {y} {z}")
                # else:
                #     sites.append(f"0 {x} {y} {z}")

    return sites


if __name__ == "__main__":
    """
    Generate standard percolation lattice data sets for benchmarking validation:
        - At theoretical pc = 0.3116 for sizes from 20 to 50 with step 5
        - From 0.2 to 0.4 with step 0.002 for sizes from 20 to 50 with step 5
    """

    # Generate percolation sites at theoretical pc
    j = 20
    probability = 0.3116
    for n in range(1):
        sites = generate_sites(probability, j)
        if not os.path.exists(f"./{j}"):
            os.makedirs(f"./{j}")
        with open(f"./{j}/percolation_sites_{probability:.4f}.xyz", "a") as f:
            f.write(f"{len(sites)}\n")
            f.write(f'Lattice="{j}.0 0.0 0.0 0.0 {j}.0 0.0 0.0 0.0 {j}.0"\n')
            for site in sites:
                f.write(site + "\n")
        f.close()
