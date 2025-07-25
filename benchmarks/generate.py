# This script generate percolation sites in a 100x100x100 orthohombic lattice

import numpy as np
from numba import jit
from typing import List
import os


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
    for j in range(10, 100, 5):
        # for i in range(1, 100):
        # generate the probability in the range [0.2, 0.4]
        probability = 0.3116
        # generate 100 configurations of the lattice
        for n in range(100):
            sites = generate_sites(probability, j)
            # if not os.path.exists(f"./{j}"):
            #     os.makedirs(f"./{j}")
            with open(f"./perco-{j}_{probability:.4f}.xyz", "a") as f:
                f.write(f"{len(sites)}\n")
                f.write(f'Lattice="{j}.0 0.0 0.0 0.0 {j}.0 0.0 0.0 0.0 {j}.0"\n')
                for site in sites:
                    f.write(site + "\n")
            f.close()
        # with open(f"./{j}/info.csv", "w") as f:
        #     f.write("project_name,probability\n")
        #     for i in range(1, 100):
        #         probability = (i * 0.002) + 0.20
        #         f.write(f"p{probability:.4f},{probability}\n")
        #     f.close()

