#!/usr/bin/env bash
# Checking typical physical-chemistry equilibrium problems from the command line.
set -e

# K from a measured degree of dissociation: N2O4 = 2 NO2, alpha = 0.201 at 298 K, 1 bar  (K ≈ 0.168)
thermoreact equilibrium "N2O4 -> 2 NO2" --init N2O4=1 --extent 0.201 -T 298.15

# K and ΔrG° from an equilibrium composition with generic species A, B, C, D  (K ≈ 9.6)
thermoreact equilibrium "A + B -> C + 2 D" --abstract --init A=2,B=1,D=3 --extent 0.79 -T 298

# van 't Hoff: K at 2000 K from K = 0.414 at 1600 K with ΔrH° = +112 kJ/mol  (K ≈ 2.23)
thermoreact vanthoff --K1 0.414 --T1 1600 --T2 2000 --dH 112

# Temperature at which K = 1, from ΔrG°(1120 K) = +22 kJ/mol and ΔrH° = +125 kJ/mol  (≈ 1359 K)
thermoreact vanthoff --dG1 22 --T1 1120 --K2 1 --dH 125

# Dehydration of copper sulfate pentahydrate: T where K = 1, from tabulated data  (≈ 398 K)
thermoreact reaction "CuSO4·5H2O -> CuSO4 + 5 H2O(g)" --crossover --cp none

# Equilibrium constant from a standard cell potential  (K ≈ 1.7e16)
thermoreact electrochem --E0 0.48 --n 2
