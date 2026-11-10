"""Physical constants (CODATA 2018 exact values) and reference conditions."""

#: Molar gas constant, J/(mol K)
R = 8.314462618
#: Faraday constant, C/mol
F = 96485.33212
#: Reference temperature for tabulated standard data, K
T_REF = 298.15
#: Standard pressure, Pa (IUPAC: 1 bar)
P_STD = 1.0e5
#: One standard atmosphere, Pa
ATM = 101325.0

#: Phase labels understood by the tool
PHASES = ("g", "l", "s", "aq")
PHASE_NAMES = {"g": "gas", "l": "liquid", "s": "solid", "aq": "aqueous (1 mol/kg)"}

#: Accepted spellings of phase/state labels -> canonical phase
PHASE_ALIASES = {
    "g": "g", "gas": "g", "v": "g", "vap": "g",
    "l": "l", "liq": "l", "liquid": "l",
    "s": "s", "cr": "s", "c": "s", "solid": "s",
    "aq": "aq", "aqueous": "aq", "ao": "aq",
}

# fmt: off
#: Standard atomic weights (IUPAC conventional values), g/mol
ATOMIC_MASSES = {
    "H": 1.008, "D": 2.014, "T": 3.016, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990,
    "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45,
    "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942,
    "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546,
    "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904,
    "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906,
    "Mo": 95.95, "Tc": 98.0, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87,
    "Cd": 112.41, "In": 114.82, "Sn": 118.71, "Sb": 121.76, "Te": 127.60, "I": 126.90,
    "Xe": 131.29, "Cs": 132.91, "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Pr": 140.91,
    "Nd": 144.24, "Pm": 145.0, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93,
    "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97,
    "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22,
    "Pt": 195.08, "Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98,
    "Po": 209.0, "At": 210.0, "Rn": 222.0, "Fr": 223.0, "Ra": 226.0, "Ac": 227.0,
    "Th": 232.04, "Pa": 231.04, "U": 238.03, "Np": 237.0, "Pu": 244.0, "Am": 243.0,
}
# fmt: on

#: Reference state of each element: (species id in the database, atoms per formula unit).
#: Used to compute formation properties at temperatures other than 298.15 K.
ELEMENT_REFERENCE = {
    "H": ("H2(g)", 2), "O": ("O2(g)", 2), "N": ("N2(g)", 2), "F": ("F2(g)", 2),
    "Cl": ("Cl2(g)", 2), "Br": ("Br2(l)", 2), "I": ("I2(s)", 2), "C": ("C(graphite)", 1),
    "S": ("S(rhombic)", 1), "P": ("P(white)", 1), "Si": ("Si(s)", 1), "He": ("He(g)", 1),
    "Ne": ("Ne(g)", 1), "Ar": ("Ar(g)", 1), "Na": ("Na(s)", 1), "K": ("K(s)", 1),
    "Mg": ("Mg(s)", 1), "Ca": ("Ca(s)", 1), "Al": ("Al(s)", 1), "Fe": ("Fe(s)", 1),
    "Cu": ("Cu(s)", 1), "Zn": ("Zn(s)", 1), "Ag": ("Ag(s)", 1), "Sn": ("Sn(s)", 1),
    "Pb": ("Pb(s)", 1), "Hg": ("Hg(l)", 1),
}
