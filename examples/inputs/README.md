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

### SiO2 — Silicon dioxide

| File | Atoms | Condition | Source / Notes |
|------|------:|-----------|----------------|
| `example-SiO2-1008at.xyz` | 1008 | <!-- e.g. ambient pressure, T = 300 K --> | <!-- origin / reference --> |
| `example-SiO2-1008at-high_pressure.xyz` | 1008 | <!-- P = ? GPa --> | |
| `example-SiO2-27216at-low_pressure.xyz` | 27216 | <!-- P = ? GPa --> | |
| `example-SiO2-27216at-11GPa.xyz` | 27216 | P = 11 GPa | |
| `example-SiO2-96000at.xyz` | 96000 | <!-- condition --> | |

**Connectivity:** `Si–O–Si` bond criterion (`"bond"` strategy)  
**Cutoffs:** Si–Si 3.50 Å · Si–O 2.30 Å · O–O 3.05 Å  
**Coordination range:** Si coordination number Z ∈ [4, 6]  
**Quickstart script:** `quickstart-a_SiO2.py`

<!-- Add any additional notes about the SiO2 dataset here. -->

---

### a-H₂O — Amorphous / high-density water

| File | Molecules | Pressure | Source / Notes |
|------|----------:|----------|----------------|
| `example-a_H2O-N8192-0.01kbar.xyz` | 8192 | 0.01 kbar | <!-- LD phase --> |
| `example-a_H2O-N8192-0.5kbar.xyz`  | 8192 | 0.5 kbar  | <!-- HD phase --> |
| `example-a_H2O-N8192-1kbar.xyz`    | 8192 | 1 kbar    | <!-- VHD phase --> |

**Connectivity:** `O–O` distance criterion (`"distance"` strategy)  
**Cutoff:** O–O 3.5 Å  
**Coordination ranges:** LD: Z = 4 · HD: Z ∈ [5, 7] · VHD: Z ≥ 8  
**Quickstart script:** `quickstart-a_H2O.py`

<!-- Add any additional notes about the a-H2O dataset here. -->

---

### a-Si — Amorphous silicon (Deringer et al 2021)

| File | Atoms | Condition | Source / Notes |
|------|------:|-----------|----------------|
| `example-a_Si-100000at-0GPa.xyz`  | 100000 | 0 GPa  | <!-- origin / reference --> |
| `example-a_Si-100000at-12GPa.xyz` | 100000 | 12 GPa | |

**Connectivity:** <!-- e.g. Si–Si distance criterion -->  
**Cutoff:** <!-- e.g. Si–Si ? Å -->  
**Quickstart script:** <!-- e.g. `quickstart-a_Si.py` (to be added) -->

<!-- Add any additional notes about the a-Si dataset here. -->

---

## Quickstart scripts

| Script | Dataset | Description |
|--------|---------|-------------|
| `quickstart-a_SiO2.py` | SiO2 | Runs pairwise coordination analysis on SiO4/SiO5/SiO6 clusters |
| `quickstart-a_H2O.py`  | a-H2O | Runs LD / HD / VHD connectivity analysis on amorphous water |

Run from the repository root:

```bash
python examples/inputs/quickstart-a_SiO2.py
python examples/inputs/quickstart-a_H2O.py
```

Output is written to `examples/outputs/`.

---

## XYZ format

All trajectory files follow the extended XYZ format:

```
<N>
Lattice="ax 0.0 0.0 0.0 ay 0.0 0.0 0.0 az" [Properties=...] [key=value ...]
<type>  <x>  <y>  <z>  [extra columns]
...
```

The `Lattice` key encodes the simulation box vectors (Å). Periodic boundary conditions are applied by default (`apply_pbc=True`).
