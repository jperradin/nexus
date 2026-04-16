# Example Input Data

This directory contains example trajectory files (XYZ format) and quickstart scripts used to demonstrate `nexus-cat` capabilities.

---

## File naming convention

```
example-<system>-<N>at[-<condition>].xyz
```

| Field         | Description                                              |
|---------------|----------------------------------------------------------|
| `<system>`    | Chemical system (e.g. `SiO2`, `a_H2O`, `a_Si`)          |
| `<N>at`       | Number of atoms (or molecules for `a_H2O`: `N<M>`)       |
| `<condition>` | Thermodynamic condition (pressure, phase label, etc.)    |

---

## Datasets

### a-SiO2 — Amorphous silica

| File | Atoms | Condition | Source / Notes |
|------|------:|-----------|----------------|
| `example-a_SiO2-8064at-0GPa.xyz` | 8064 | T = 300K, P = 0GPa | Perradin et al. PhD thesis 2025 |
| `example-a_SiO2-8064at-5GPa.xyz` | 8064 | T = 300K, P = 5GPa | |
| `example-a_SiO2-8064at-10GPa.xyz` | 8064 | T = 300K, P = 10GPa | |
| `example-a_SiO2-8064at-15GPa.xyz` | 8064 | T = 300K, P = 15GPa | |
| `example-a_SiO2-8064at-20GPa.xyz` | 8064 | T = 300K, P = 20GPa | |
| `example-a_SiO2-8064at-25GPa.xyz` | 8064 | T = 300K, P = 25GPa | |
| `example-a_SiO2-8064at-30GPa.xyz` | 8064 | T = 300K, P = 30GPa | |
| `example-a_SiO2-8064at-35GPa.xyz` | 8064 | T = 300K, P = 35GPa | |

**Connectivity:** `Si–O–Si` (`"CoordinationStrategy"` + `"SharedStrategy`)  
**Cutoffs:** Si–Si 3.50 Å · Si–O 2.30 Å · O–O 3.05 Å  
**Coordination range:** Si coordination number Z ∈ [4, 6]
**Quickstart script:** `quickstart-a_SiO2.py`

  1) First run :
    - SiO4-SiO4
    - SiO4-SiO5
    - SiO5-SiO5
    - SiO5-SiO6
    - SiO6-SiO6
  2) Second run :
    - Stishovite (SiO6*-SiO6)

---

### a-H₂O — Amorphous ice

| File | Molecules | Pressure | Source / Notes |
|------|----------:|----------|----------------|
| `example-a_H2O-8192mol-0.01kbar.xyz` | 8192 | T = 124K, P = 0.01 kbar | Hasmy et al. 2025 |
| `example-a_H2O-8192mol-0.5kbar.xyz`  | 8192 | T = 124K, P = 0.5 kbar  | - |
| `example-a_H2O-8192mol-1kbar.xyz`    | 8192 | T = 124K, P = 1 kbar    | - |

**Connectivity:** `O–O` distance criterion (`"distance"` strategy)  
**Cutoff:** O–O 3.5 Å  
**Coordination ranges:** LD: Z = 4 · HD: Z ∈ [5, 7] · VHD: Z ≥ 8  
**Quickstart script:** `quickstart-a_H2O.py`

  1) First run : LD
  2) Second run : HD
  3) Third run : VHD

---

### a-Si — Amorphous silicon (Deringer et al 2021)

| File | Atoms | Condition | Source / Notes |
|------|------:|-----------|----------------|
| `example-a_Si-100000at-0GPa.xyz`  | 100000 | 0 GPa  | Deringer et al. 2021 |
| `example-a_Si-100000at-12GPa.xyz` | 100000 | 12 GPa | |

**Connectivity:** `Si–Si` distance criterion (`"distance"` strategy)  
**Cutoff:** Si–Si 3.0 Å  
**Coordination ranges:** LD: Z = 4 · HD: Z ∈ [5, 7] · VHD: Z ≥ 8  
**Quickstart script:** `quickstart-a_Si.py`

  1) First run : LD
  2) Second run : HD
  3) Third run : VHD

---

## Quickstart scripts

| Script | Dataset | Description |
|--------|---------|-------------|
| `quickstart-a_SiO2.py` | SiO2 | Runs pairwise coordination analysis on SiO4/SiO5/SiO6 clusters |
| `quickstart-a_H2O.py`  | a-H2O | Runs LD / HD / VHD connectivity analysis on amorphous water |
| `quickstart-a_Si.py`  | a-Si | Runs LD / HD / VHD connectivity analysis on amorphous silicon |

Run from the examples inputs directory:

```bash
cd examples/inputs
python quickstart-a_SiO2.py # or
python quickstart-a_H2O.py # or
python quickstart-a_Si.py
```

Output is written to `examples/outputs/`.

---

## XYZ format

All trajectory files follow the extended XYZ format:

```
<N>
Lattice="lx 0.0 0.0 0.0 ly 0.0 0.0 0.0 lz" 
<type>  <x>  <y>  <z>  
...
```

The `Lattice` key encodes the simulation box vectors (Å). 
Periodic boundary conditions are applied by default (`apply_pbc=True`).

