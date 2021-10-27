import argparse
from typing import List, Dict, Tuple, TextIO
from structures import Atom, Molecule


def skip_two_lines(file: TextIO) -> None:
    for _ in range(2):
        file.readline()


def atoms_bonds_count(file: TextIO) -> (int, int):
    """Returns number of atoms and number of bonds within a molecule"""
    skip_two_lines(file)
    line = file.readline()
    # AAABBB
    return int(line[:3]), int(line[3:6])


def get_atoms(file: TextIO, number_of_atoms: int) -> List[Atom]:
    """Saves data about atoms in sdf file."""
    atoms = []
    for _ in range(number_of_atoms):
        x, y, z, element_name, *rest = file.readline().split()
        x, y, z = float(x), float(y), float(z)

        atom = Atom(x, y, z, element_name)
        atoms.append(atom)

    return atoms


def bonds(file: TextIO, number_of_bonds: int) -> Dict[Tuple[int, int], int]:
    """Process bonds part of sdf file."""
    bonds_data = {}
    for i in range(number_of_bonds):
        line = file.readline()
        first_atom, second_atom = sorted((int(line[:3]), int(line[3:6])))
        type_of_bond = int(line[6:9])

        bonds_data[(first_atom, second_atom)] = type_of_bond

    return bonds_data


def skip_non_structural_data(file: TextIO, delimiter: str, pdb: bool = False) -> str:
    """Skip lines until the line with delimiter is reached.
    If pdb=True - usually searching for delimiter MODEL.
    But it is optional - in case that there is no MODEL in pdb,
    program returns line with keyword ATOM"""
    line = file.readline()
    while not line.startswith(delimiter):
        if not pdb:
            line = file.readline()
            continue
        else:
            if line.startswith('ATOM'):
                break
            line = file.readline()
    return line


def find_min_max(molecules: List[Molecule]) -> (Molecule, Molecule):
    """Returns largest and smallest molecules."""
    max_molecule = molecules[0]
    min_molecule = molecules[0]
    for molecule in molecules:
        if len(molecule) > len(max_molecule):
            max_molecule = molecule
        elif len(molecule) < len(min_molecule):
            min_molecule = molecule
    return max_molecule, min_molecule


def statistics(molecules: List[Molecule]) -> None:
    """Prints out number of all molecules, the largest molecule and
        smallest molecule.
    """
    max_molecule, min_molecule = find_min_max(molecules)

    print(f'Number of molecules: {len(molecules)}')
    print(f'Highest number of atoms has molecule: '
          f'{max_molecule.name}'
          f', exactly: {len(max_molecule.atoms)} atoms.')
    print(f'Smallest number of atoms has molecule: '
          f'{min_molecule.name}'
          f', exactly: {len(min_molecule.atoms)} atoms.')


def read_sdf_file(path_to_file: str) -> None:
    with open(path_to_file) as file:
        molecules = []
        molecular_name = file.readline().strip()
        while molecular_name:
            number_of_atoms, number_of_bonds = atoms_bonds_count(file)
            atoms = get_atoms(file, number_of_atoms)
            bonds_data = bonds(file, number_of_bonds)
            skip_non_structural_data(file, '$$$$\n')

            molecule = Molecule(molecular_name, atoms, bonds_data)

            # print(molecule.name)
            # molecule.molecular_formula()
            # print()

            molecules.append(molecule)

            molecular_name = file.readline().strip()

        statistics(molecules)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='Input file in .sdf format')
    args = parser.parse_args()

    input_file = args.file
    read_sdf_file(input_file)
