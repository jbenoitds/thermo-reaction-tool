"""thermoreact - reaction thermochemistry and equilibrium calculator.

Quick start::

    >>> from thermoreact import Reaction
    >>> rxn = Reaction.from_string("N2 + 3 H2 -> 2 NH3")
    >>> p = rxn.properties(298.15)
    >>> round(p.dH / 1000, 1), round(p.dG / 1000, 1)
    (-91.8, -32.7)
"""
from .constants import F, P_STD, R, T_REF
from .database import AmbiguousSpecies, Database, SpeciesNotFound
from .equilibrium import EquilibriumResult, k_from_extent, solve_equilibrium
from .formula import FormulaError, hill_formula, molar_mass, parse_formula
from .joback import JobackError, estimate
from .reaction import Reaction, ReactionError, ReactionProperties
from .smiles import Molecule, SmilesError, parse_smiles
from .species import Species

__version__ = "0.1.0"

__all__ = [
    "R", "F", "P_STD", "T_REF",
    "Database", "Species", "SpeciesNotFound", "AmbiguousSpecies",
    "Reaction", "ReactionProperties", "ReactionError",
    "solve_equilibrium", "k_from_extent", "EquilibriumResult",
    "parse_formula", "hill_formula", "molar_mass", "FormulaError",
    "parse_smiles", "Molecule", "SmilesError", "estimate", "JobackError",
    "__version__",
]
