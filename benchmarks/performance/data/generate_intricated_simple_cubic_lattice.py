import numpy as np
import os

with open("dummy_sites.xyz", "w") as f:
    f.write("54000\n")
    f.write('Lattice="30.0 0.0 0.0 0.0 30.0 0.0 0.0 0.0 30.0"\n')
    for i in range(30):
        for j in range(30):
            for k in range(30):
                f.write(f"1 {i} {j} {k}\n")
                f.write(f"2 {i+0.5} {j+0.5} {k+0.5}\n")

f.close()
