import numpy as np
import os
from tqdm import tqdm
from numba import njit, prange


# @njit(parallel=True)
def generate_sites(num_sites, lattice_size):
    networking_sites = []
    bridge_sites = []
    network_hash = set()
    bridge_hash = set()
    limit = np.int32(lattice_size**3)

    for _ in tqdm(range(num_sites)):
        choice = np.random.randint(0, 2)
        if choice == 0 and len(bridge_sites) < limit:
            while True:
                x = np.random.randint(0, lattice_size)
                y = np.random.randint(0, lattice_size)
                z = np.random.randint(0, lattice_size)
                pos = np.array([x + 0.5, y + 0.5, z + 0.5], dtype=np.float32)
                hash = f"{x}{y}{z}"
                if hash not in bridge_hash:
                    bridge_hash.add(hash)
                    bridge_sites.append(pos)
                    break
        else:
            while True:
                x = np.random.randint(0, lattice_size)
                y = np.random.randint(0, lattice_size)
                z = np.random.randint(0, lattice_size)
                pos = np.array([x, y, z], dtype=np.float32)
                hash = f"{x}{y}{z}"
                if hash not in network_hash:
                    network_hash.add(hash)
                    networking_sites.append(pos)
                    break

    return networking_sites, bridge_sites


if __name__ == "__main__":
    """
    Generate simple cubic lattice data as two intricated percolation lattice 1 and 2.
    Sets generated for benchmarking performance:
     N - Lattice size - concentration
     50000 - 50 - 0.4
     75000 - 50 - 0.6
     86400 - 60 - 0.4
     129600 - 60 - 0.6
     137200 - 70 - 0.4
     205800 - 70 - 0.6
     204800 - 80 - 0.4
     307200 - 80 - 0.4
     291600 - 90 - 0.6
     437400 - 90 - 0.6
     400000 - 100 - 0.4
     600000 - 100 - 0.6
    """
    lattices = [
        (50000, 50),
        (75000, 50),
        (86400, 60),
        (129600, 60),
        (137200, 70),
        (205800, 70),
        (204800, 80),
        (307200, 80),
        (291600, 90),
        (437400, 90),
        (400000, 100),
        (600000, 100),
    ]

    lattice_size = 30  # Size of the cubic lattice
    num_sites = 16200  # Total number of sites to generate, N < lattice_size^3

    # for num_sites, lattice_size in lattices:

    networking_sites, bridge_sites = generate_sites(num_sites, lattice_size)

    n_sites = len(networking_sites) + len(bridge_sites)
    lattice_string = f'Lattice="{lattice_size} 0.0 0.0 0.0 {lattice_size} 0.0 0.0 0.0 {lattice_size}"'

    with open(f"./benchmark-{num_sites}-{lattice_size}.xyz", "w") as f:
        f.write(f"{n_sites}\n")
        f.write(f"{lattice_string}\n")
        for site in networking_sites:
            f.write(f"1 {site[0]} {site[1]} {site[2]}\n")
        for site in bridge_sites:
            f.write(f"2 {site[0]} {site[1]} {site[2]}\n")
