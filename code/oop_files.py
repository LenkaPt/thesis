import argparse
from typing import List, Dict, Tuple, TextIO
from collections import Counter


class Atom:
    def __init__(self, x, y, z, name):
        self._x = x
        self._y = y
        self._z = z
        self._name = name

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    @property
    def name(self):
        return self._name

    def print_coordinates(self):
        print(f'Coordinates of atom {self._name}: \n\t'
              f'x: {self._x}, y: {self._y}, z: {self._z}')


class Molecule:
    def __init__(self, name, atoms, bonds):
        self._name = name
        self._atoms = atoms
        self._bonds = bonds

    @property
    def name(self):
        return self._name

    @property
    def atoms(self):
        return self._atoms

    def molecular_formula(self):
        atoms = [atom.name for atom in self._atoms]
        for atom, count in Counter(atoms).items():
            print(f'{atom}{count}', end='')
        print()

    def __len__(self):
        return len(self.atoms)


'''atom1 = Atom(1, 2, 3, 'H')
atom2 = Atom(0.1, 0.2, 0.3, 'H')
atom3 = Atom(10, 20, 30, 'O')
atom4 = Atom(1, 2, 3, 'O')

molec = Molecule('H2O', [atom1, atom2, atom3, atom4], {(1, 12): 1, (1, 17): 1, (1, 18): 1, (2, 16): 2, (3, 17): 2})
molec.molecular_formula()
print(molec.atoms)'''


class Chain:
    def __init__(self, residues, name):
        self.residues = residues
        self._name = name


class Residue:
    def __init__(self, name, number, atoms):
        self._name = name
        self._number = number
        self._atoms = atoms

    @property
    def name(self):
        return self._name


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


def skip_non_structural_data(file: TextIO, delimiter: str) -> str:
    """Skip lines until the line with delimiter is reached."""
    line = file.readline()
    while not line.startswith(delimiter):
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

            '''print(molecule.name)
            molecule.molecular_formula()
            print()'''

            molecules.append(molecule)

            molecular_name = file.readline().strip()

        statistics(molecules)


######################
def get_atom_pdb(line: str) -> Atom:
    """Returns one atom from pdb file (as instance of Atom class)"""
    atom_name = line[13:16]
    x, y, z = float(line[31:38]), float(line[39:46]), float(line[47:54])
    atom = Atom(x, y, z, atom_name)
    return atom


def get_residue_pdb(file: TextIO, line: str) -> Tuple[Residue, str, str]:
    """Returns one residuum from pdb file."""
    chain_id = line[22]
    residue_seq_number = line[23:26]
    residue_name = line[17:20]
    atoms = []
    while residue_seq_number == line[23:26]:
        atom = get_atom_pdb(line)
        atoms.append(atom)
        line = file.readline()
    residue = Residue(residue_name, residue_seq_number, atoms)

    return residue, chain_id, line


def read_pdb_file(path_to_file: str) -> None:
    with open(path_to_file) as file:
        line = skip_non_structural_data(file, 'ATOM')
        residues = []
        for i in range(20):
            residue, chain_id, line = get_residue_pdb(file, line)
            residues.append(residue)
            print(residues[i].name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='Input file in .sdf or .pdb format')
    args = parser.parse_args()

    input_file = args.file
    read_pdb_file(input_file)
