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
    # AAABBB (integers)
    try:
        return int(line[:3]), int(line[3:6])
    except ValueError:
        raise ValueError('Please, check the counts line (fourth line). '
                         'Expected format is aaabbb, where aaa is '
                         'number of atoms (integer) '
                         'and bbb is number of bonds (integer).')


def get_atoms(file: TextIO, number_of_atoms: int) -> List[Atom]:
    """Saves data about atoms in sdf file."""
    atoms = []
    for _ in range(number_of_atoms):
        line = file.readline()
        x, y, z, element_name = line[:10], line[10:20], line[20:30], line[30:34].strip()
        # x, y, z should be floats
        try:
            x, y, z = float(x), float(y), float(z)
        except ValueError:
            if not element_name:
                raise ValueError(f'Please check if number of atoms in '
                                 f'the counts line is correct. '
                                 f'Coordinate information about atom are expected.'
                                 f'\nGiven line: {line.strip()}')
            raise ValueError(f'Please check coordinates of atom number '
                             f'{len(atoms) + 1} ({element_name}). '
                             f'Coordinates must be floats.')

        atom = Atom(x, y, z, element_name)
        atoms.append(atom)

    return atoms


def bonds(file: TextIO, number_of_bonds: int) -> Dict[Tuple[int, int], int]:
    """Process bonds part of sdf file."""
    bonds_data = {}
    for i in range(number_of_bonds):
        line = file.readline()
        try:
            first_atom, second_atom = sorted((int(line[:3]), int(line[3:6])))
            type_of_bond = int(line[6:9])
        except ValueError:
            if '.' in line[3:6]:
                # probably info about atom coordinates expected
                raise ValueError(f'Please check if number of atoms in '
                                 f'the counts line is correct.'
                                 f' Bond block expected.'
                                 f'\nGiven line: {line.strip()}')
            raise ValueError(f'Please check {len(bonds_data) + 1}. line '
                             f'of the bond block. Format 111222ttt is '
                             f'expected, where all numbers must be integers. '
                             f'111 is number of first atom, 222 is number of second'
                             f' atom and ttt is bond type.')

        bonds_data[(first_atom, second_atom)] = type_of_bond

    return bonds_data


def skip_non_structural_data(file: TextIO, delimiters: List[str]) -> str:
    """Skip lines until the line with one of the delimiters is reached."""
    delimiters = tuple(delimiters)
    line = file.readline()
    while not line.startswith(delimiters):
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
            try:
                number_of_atoms, number_of_bonds = atoms_bonds_count(file)
                atoms = get_atoms(file, number_of_atoms)
                bonds_data = bonds(file, number_of_bonds)
            except ValueError as e:
                raise ValueError(f'Problem occurred! \nMolecule: {molecular_name}'
                                 f'\n{e}\n-------')
            skip_non_structural_data(file, ['$$$$\n'])

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
