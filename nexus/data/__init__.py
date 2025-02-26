import numpy as np

__all__ = [
    "correlation_lengths",
    "atomic_masses",
    "atomic_numbers",
    "chemical_symbols",
]

# Correlation lengths for elements (in fm)
correlation_lengths = np.array(
    [
        -3.7390,  # H
        np.nan,  # He
        -1.90,  # Li
        7.79,  # Be
        5.30,  # B
        6.6460,  # C
        9.36,  # N
        5.803,  # O
        5.654,  # F
        4.566,  # Ne
        3.63,  # Na
        5.375,  # Mg
        3.449,  # Al
        4.149,  # Si
        5.13,  # P
        2.847,  # S
        9.5770,  # Cl
        1.909,  # Ar
        3.67,  # K
        4.70,  # Ca
        12.29,  # Sc
        -3.438,  # Ti
        -0.3824,  # V
        3.635,  # Cr
        -3.73,  # Mn
        9.45,  # Fe
        2.49,  # Co
        10.3,  # Ni
        7.718,  # Cu
        5.680,  # Zn
        7.288,  # Ga
        8.185,  # Ge
        6.58,  # As
        7.970,  # Se
        6.795,  # Br
        7.81,  # Kr
        7.09,  # Rb
        7.02,  # Sr
        7.75,  # Y
        7.16,  # Zr
        7.054,  # Nb
        6.715,  # Mo
        np.nan,  # Tc
        7.03,  # Ru
        5.88,  # Rh
        5.91,  # Pd
        5.922,  # Ag
        4.87,  # Cd
        4.065,  # In
        6.225,  # Sn
        5.57,  # Sb
        5.80,  # Te
        5.28,  # I
        4.92,  # Xe
        5.42,  # Cs
        5.07,  # Ba
        8.24,  # La
        4.84,  # Ce
        4.58,  # Pr
        7.69,  # Nd
        np.nan,  # Pm
        0.80,  # Sm
        7.22,  # Eu
        6.5,  # Gd
        7.38,  # Tb
        16.9,  # Dy
        8.01,  # Ho
        7.79,  # Er
        7.07,  # Tm
        12.43,  # Yb
        7.21,  # Lu
        7.7,  # Hf
        6.91,  # Ta
        4.86,  # W
        9.2,  # Re
        10.7,  # Os
        10.6,  # Ir
        9.60,  # Pt
        7.63,  # Au
        12.692,  # Hg
        12.692,  # Tl
        9.405,  # Pb
        8.532,  # Bi
        np.nan,  # Po
        np.nan,  # At
        np.nan,  # Rn
        np.nan,  # Fr
        np.nan,  # Ra
        np.nan,  # Ac
        10.31,  # Th
        9.1,  # Pa
        8.417,  # U
        np.nan,  # Np
        np.nan,  # Pu
        np.nan,  # Am
        np.nan,  # Cm
        np.nan,  # Bk
        np.nan,  # Cf
        np.nan,  # Es
        np.nan,  # Fm
        np.nan,  # Md
        np.nan,  # No
        np.nan,  # Lr
        np.nan,  # Rf
        np.nan,  # Db
        np.nan,  # Sg
        np.nan,  # Bh
        np.nan,  # Hs
        np.nan,  # Mt
        np.nan,  # Ds
        np.nan,  # Rg
        np.nan,  # Cn
        np.nan,  # Nh
        np.nan,  # Fl
        np.nan,  # Mc
        np.nan,  # Lv
        np.nan,  # Ts
        np.nan,  # Og
    ]
)

# Chemical symbols for elements
chemical_symbols = np.array(
    [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
        "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
        "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
        "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc",
        "Lv", "Ts", "Og",
    ]
)

# Atomic numbers for elements
atomic_numbers = {symbol: Z for Z, symbol in enumerate(chemical_symbols)}

# Atomic masses for elements (IUPAC 2016)
atomic_masses = np.array(
    [
        1.008, 4.002602, 6.94, 9.0121831, 10.81, 12.011, 14.007, 15.999, 18.998403163, 20.1797, 22.98976928, 24.305,
        26.9815385, 28.085, 30.973761998, 32.06, 35.45, 39.948, 39.0983, 40.078, 44.955908, 47.867, 50.9415, 51.9961,
        54.938044, 55.845, 58.933194, 58.6934, 63.546, 65.38, 69.723, 72.630, 74.921595, 78.971, 79.904, 83.798, 85.4678,
        87.62, 88.90584, 91.224, 92.90637, 95.95, 97.90721, 101.07, 102.90550, 106.42, 107.8682, 112.414, 114.818, 118.710,
        121.760, 127.60, 126.90447, 131.293, 132.90545196, 137.327, 138.90547, 140.116, 140.90766, 144.242, 144.91276,
        150.36, 151.964, 157.25, 158.92535, 162.500, 164.93033, 167.259, 168.93422, 173.054, 174.9668, 178.49, 180.94788,
        183.84, 186.207, 190.23, 192.217, 195.084, 196.966569, 200.592, 204.38, 207.2, 208.98040, 208.98243, 209.98715,
        222.01758, 223.01974, 226.02541, 227.02775, 232.0377, 231.03588, 238.02891, 237.04817, 244.06421, 243.06138,
        247.07035, 247.07031, 251.07959, 252.0830, 257.09511, 258.09843, 259.1010, 262.110, 267.122, 268.126, 271.134,
        270.133, 269.1338, 278.156, 281.165, 281.166, 285.177, 286.182, 289.190, 289.194, 293.204, 293.208, 294.214,
    ]
)
