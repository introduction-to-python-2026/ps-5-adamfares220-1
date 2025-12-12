# Add the import statements for functions from string_utils.py and equation_utils.py here
from string_utils import parse_chemical_reaction, count_atoms_in_reaction
from equation_utils import build_equations, my_solve


import re


def count_atoms(compound):
    return {el: int(count) if count else 1 
            for el, count in re.findall(r'([A-Z][a-z]?)([0-9]*)', compound)}

def parse_reaction_simple(reaction_string):
    reactants_str, products_str = reaction_string.replace(" ", "").split('->')
    reactants = [count_atoms(c) for c in reactants_str.split('+')]
    products = [count_atoms(c) for c in products_str.split('+')]
    return reactants, products


def balance_reaction(reaction): 
 
    reactant_atoms, product_atoms = parse_reaction_simple(reaction)
    
    print( {reactant_atoms})
    print({product_atoms})
