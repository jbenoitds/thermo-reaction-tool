"""Using thermoreact from Python instead of the command line.

Run with:  python examples/python_api.py
"""
from thermoreact import Database, Reaction, solve_equilibrium

# 1. Reaction properties at several temperatures --------------------------------------
rxn = Reaction.from_string("CO + H2O(g) -> CO2 + H2")          # water-gas shift
for T in (298.15, 600, 1000):
    p = rxn.properties(T)
    print("T = %6.1f K  ΔrH° = %7.2f kJ/mol  ΔrG° = %7.2f kJ/mol  K = %.4g" % (T, p.dH / 1e3, p.dG / 1e3, p.K))

# 2. Automatic balancing and a crossover temperature -----------------------------------
rxn = Reaction.from_string("CaCO3 -> CaO + CO2", balance=True)
print("\n%s : K = 1 at %.0f K" % (rxn, rxn.crossover_temperature()[0]))

# 3. A molecule that is not in the database: given as SMILES, estimated (Joback) -------
rxn = Reaction.from_string("smiles:CC(C)CC + O2 -> CO2 + H2O(l)", balance=True)
print("\n%s\n  ΔcH° = %.0f kJ/mol (estimated)" % (rxn, rxn.properties().dH / 1e3))

# 4. Equilibrium composition: ammonia synthesis at 700 K and 200 bar -------------------
rxn = Reaction.from_string("N2 + 3 H2 -> 2 NH3")
eq = solve_equilibrium(rxn, {"N2": 1, "H2": 3}, T=700, P_bar=200)
print("\nNH3 mole fraction at 700 K, 200 bar: %.3f" % eq.mole_fractions()["NH3(g)"])

# 5. Species data -----------------------------------------------------------------------
db = Database.default()
water = db.get("H2O(g)")
print("\nH2O(g): Cp(1000 K) = %.2f J/(mol K), ΔfG°(1000 K) = %.2f kJ/mol"
      % (water.cp(1000), db.formation(water, 1000)["dGf"] / 1e3))
