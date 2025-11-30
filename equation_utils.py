import re
from collections import defaultdict
import numpy as np
from sympy import Matrix, lcm 

def parse_chemical_reaction(reaction_string):
    sides = reaction_string.replace(" ", "").split('->')
    reactants = sides.split('+')
    products = sides.split('+')
    return reactants, products

def count_atoms_in_reaction(compounds):
    atom_counts_list = []
    for compound in compounds:
        counts = defaultdict(int)
        for element, count in re.findall(r'([A-Z][a-z]?)([0-9]*)', compound):
            counts[element] += int(count) if count else 1
        atom_counts_list.append(counts)
    return atom_counts_list

def build_equations(reactant_atoms, product_atoms):
    all_elements = sorted(set().union(*reactant_atoms, *product_atoms))
    num_compounds = len(reactant_atoms) + len(product_atoms)
    num_elements = len(all_elements)
    
    matrix = np.zeros((num_elements, num_compounds), dtype=int)
    
    for j, counts in enumerate(reactant_atoms):
        for i, element in enumerate(all_elements):
            matrix[i, j] = counts[element]
            
    for j, counts in enumerate(product_atoms):
        for i, element in enumerate(all_elements):
            matrix[i, len(reactant_atoms) + j] = -counts[element]
            
    return matrix

def my_solve(equations_matrix):
    A = Matrix(equations_matrix)
    coeffs = A.nullspace()
    
    denominator_lcm = lcm([term.q for term in coeffs])
    integer_coeffs = [int(term * denominator_lcm) for term in coeffs]
    
    scaling_factor = coeffs / integer_coeffs
    
    return [scaling_factor * int(term) for term in integer_coeffs]

def balance_reaction(reaction):
    reactants_str, products_str = parse_chemical_reaction(reaction)
    reactant_atoms = count_atoms_in_reaction(reactants_str)
    product_atoms = count_atoms_in_reaction(products_str)

    equations_matrix = build_equations(reactant_atoms, product_atoms)
    coefficients = my_solve(equations_matrix)

    return coefficients


ELEMENTS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'I', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo'
]

def generate_equation_for_element(compounds, coefficients, element):
    """Generates a symbolic equation for the given element from compounds and coefficients.  
    Example: For H in reactants [{'H': 2}, {'O': 4, 'H': 1}], coefficients [a0, a1], returns 2*a0 + a1."""
    equation = 0
    for i, compound in enumerate(compounds):
        if element in compound:
            equation += coefficients[i] * compound[element]
    return equation


def build_equations(reactant_atoms, product_atoms):
    """Builds a list of symbolic equations for each element to balance a chemical reaction.  
    Example: For H2 + O2 -> H2O, returns equations [2*a0 - 2*b0, a1 - b0]."""
    ## coefficients ##
    reactant_coefficients = list(symbols(f'a0:{len(reactant_atoms)}'))
    product_coefficients = list(symbols(f'b0:{len(product_atoms)}')) 
    product_coefficients = product_coefficients[:-1] + [1] # Ensure the last coefficient is 1

    ## equations ##
    equations = []
    for element in ELEMENTS:
        lhs = generate_equation_for_element(reactant_atoms, reactant_coefficients, element)
        rhs = generate_equation_for_element(product_atoms, product_coefficients, element)
        if lhs != 0 or rhs != 0:
            equations.append(Eq(lhs, rhs))

    return equations, reactant_coefficients + product_coefficients[:-1]


def my_solve(equations, coefficients):
    """Solves the system of equations for the coefficients of the reaction.  
    Example: For equations [2*a0 - 2*b0, a1 - b0], returns [1.0, 1.0]."""
    solution = sympy_solve(equations, coefficients)

    if len(solution) == len(coefficients):
        coefficient_values = list()
        for coefficient in coefficients:
            coefficient_values.append(float(solution[coefficient]))
        return coefficient_values





