def split_before_each_uppercases(formula):
   if not formula:
       return []
   start = 0
   split_formula = []
   for end in range(1, len(formula)):
       if formula[end].isupper():
           split_formula.append(formula[start:end])
           start = end

   split_formula.append(formula[start:len(formula)])

  return split_formula

def split_at_first_digit(formula):
    digit_location = 0
    for char in formula[1:]:
        digit_location += 1
        if char.isdigit():
            break
    else:
        digit_location = len(formula)

    if digit_location == len(formula):
        return (formula, 1)
    else:
        prefix = formula[:digit_location]
        numeric_part_str = formula[digit_location:]
        return (prefix, int(numeric_part_str))

def count_atoms_in_molecule(molecular_formula):
    atoms_count = {}
    for atom in split_by_capitals(molecular_formula):
        atom_name, atom_count = split_at_number(atom)
        if atom_name in atoms_count:
            atoms_count[atom_name] += int(atom_count) if atom_count else 1
        else:
            atoms_count[atom_name] = int(atom_count) if atom_count else 1
    return atoms_count

def parse_chemical_reaction(reaction_equation):
    reaction_equation = reaction_equation.replace(" ", "")
    if '->' in reaction_equation:
        parts = reaction_equation.split('->')
        reactants_str = parts
        products_str = parts
    else:
        return ([], [])
    reactants = reactants_str.split('+')
    products = products_str.split('+')
    return (reactants, products)



def parse_chemical_reaction(reaction_equation):
    """Takes a reaction equation (string) and returns reactants and products as lists.  
    Example: 'H2 + O2 -> H2O' → (['H2', 'O2'], ['H2O'])"""
    reaction_equation = reaction_equation.replace(" ", "")  # Remove spaces for easier parsing
    reactants, products = reaction_equation.split("->")
    return reactants.split("+"), products.split("+")

def count_atoms_in_reaction(molecules_list):
    """Takes a list of molecular formulas and returns a list of atom count dictionaries.  
    Example: ['H2', 'O2'] → [{'H': 2}, {'O': 2}]"""
    molecules_atoms_count = []
    for molecule in molecules_list:
        molecules_atoms_count.append(count_atoms_in_molecule(molecule))
    return molecules_atoms_count
