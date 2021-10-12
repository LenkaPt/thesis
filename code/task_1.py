"""First task.
17.9.2021 - 23.9.2021
Load and process sdf file.
"""

import argparse
from typing import List, Union, Dict, Tuple, TextIO
from collections import Counter, defaultdict


def skip_two_lines(file: TextIO) -> None:
    for _ in range(2):
        file.readline()


def atoms_bonds_count(file: TextIO) -> (int, int):
    """Returns number of atoms and number of bonds within a molecule"""
    skip_two_lines(file)
    line = file.readline()
    # AAABBB
    return int(line[:3]), int(line[3:6])


def print_info_about_molecule(molecule_name: str, number_of_atoms: int,
                              number_of_bonds: int) -> None:
    """Prints name of molecule, number of its atoms
     and number of its bonds.
     """
    print(f'Molecular name: {molecule_name}')
    print(f'Number of atoms: {number_of_atoms}.\n'
          f'Number of bonds: {number_of_bonds}.\n')


def coordinates(file: TextIO,
                number_of_atoms: int) -> List[Dict[str, Union[float, str]]]:
    """Process coordinates part of file.
    Info about one atom is stored in dictionary.
    """
    coordinates_of_atoms = []
    for _ in range(number_of_atoms):
        atom = {}
        x, y, z, element_name, *rest = file.readline().split()
        x, y, z = float(x), float(y), float(z)

        atom['x'], atom['y'], atom['z'], atom['element name'] = x, y, z, element_name
        coordinates_of_atoms.append(atom)

    return coordinates_of_atoms


def bonds(file: TextIO, number_of_bonds: int) -> Dict[Tuple[int, int], int]:
    """Process bonds part of file.
    Info about one bond is stored in tuple.
    """
    bonds_data = {}
    for i in range(number_of_bonds):
        line = file.readline()
        first_atom, second_atom = sorted((int(line[:3]), int(line[3:6])))
        type_of_bond = int(line[6:9])

        bonds_data[(first_atom, second_atom)] = type_of_bond

    return bonds_data


def molecular_formula(number_of_atoms: int,
                      molecule_data: List[Dict[str, Union[float, str]]]) -> None:
    """Prints molecular formula."""
    atom_names_list = [molecule_data[i]['element name'] for i in range(number_of_atoms)]
    print('Molecular formula: ', end='')
    for key, value in Counter(atom_names_list).items():
        print(f'{key}{value}', end=' ')
    print()


def non_structural_data(file: TextIO, delimiter: str) -> str:
    """Skip lines until the line with delimiter is reached."""
    line = file.readline()
    while not line.startswith(delimiter):
        line = file.readline()
    return line


def bonds_of_atom(atoms_bonds_data: Dict[int, Counter],
                  atom: Union[int, str]) -> None:
    """Prints information about bonds outgoing of specific atom."""
    print(f'Atom {atom} forms {atoms_bonds_data[atom][1]}',
          f'single bonds, {atoms_bonds_data[atom][2]} double bonds and ',
          f'{atoms_bonds_data[atom][3]} triple bonds.')


def process_data_about_bonds(
        data_about_bonds: Dict[Tuple[Union[int, str], Union[int, str]], int]) \
        -> Dict[int, Counter]:
    """Provides informations about atoms and its bonds.
    Returns dictionary. Key -> number of atom or atom name.
    Value is also dictionary, where key is type of bond and value is
    how many bonds of that type forms particular atom.
    """
    atoms_bonds = defaultdict(Counter)

    for atoms, bond in data_about_bonds.items():
        first_atom, sec_atom = atoms
        atoms_bonds[first_atom][bond] += 1
        atoms_bonds[sec_atom][bond] += 1

    bonds_of_atom(atoms_bonds, 'CG1')
    bonds_of_atom(atoms_bonds, 8)
    bonds_of_atom(atoms_bonds, 'C')

    return atoms_bonds


def find_min_max(sequence: List[int]) -> (int, int):
    """Returns index of max and index of min element"""
    max_index = 0
    min_index = 0
    for i in range(len(sequence)):
        if sequence[i] > sequence[max_index]:
            max_index = i
        elif sequence[i] < sequence[min_index]:
            min_index = i
    return max_index, min_index


def statistics(names_of_all_molecules: List[str],
               number_of_atoms_of_each_molecule: List[int]) -> None:
    """Prints out number of all molecules, the largest molecule and
    smallest molecule.
    """
    max_index, min_index = find_min_max(number_of_atoms_of_each_molecule)

    print(f'Number of molecules: {len(number_of_atoms_of_each_molecule)}')
    print(f'Highest number of atoms has molecule: '
          f'{names_of_all_molecules[max_index]}'
          f', exactly: {number_of_atoms_of_each_molecule[max_index]} atoms.')
    print(f'Smallest number of atoms has molecule: '
          f'{names_of_all_molecules[min_index]}'
          f', exactly: {number_of_atoms_of_each_molecule[min_index]} atoms.')


def read_sdf_file(path_to_file: str) -> None:
    with open(path_to_file) as file:
        names_of_all_molecules = []
        number_of_atoms_of_each_molecule = []

        while True:
            molecule_name = file.readline().strip()
            if not molecule_name:
                break

            names_of_all_molecules.append(molecule_name)

            number_of_atoms, number_of_bonds = atoms_bonds_count(file)
            number_of_atoms_of_each_molecule.append(number_of_atoms)

            print_info_about_molecule(molecule_name, number_of_atoms,
                                      number_of_bonds)

            coordinates_data = coordinates(file, number_of_atoms)
            data_about_bonds = bonds(file, number_of_bonds)
            process_data_about_bonds(data_about_bonds)

            molecular_formula(number_of_atoms, coordinates_data)
            non_structural_data(file, '$$$$\n')

            print('--------\n')
        statistics(names_of_all_molecules, number_of_atoms_of_each_molecule)


###############################################
"""Task 2
24.9.2021 - 5.10.2021
Load and process pdb file.
"""


def atom_coordinates_pdb(line: str) -> Dict[str, Union[float, str]]:
    """Returns coordinates of atom (pdb format)."""
    atom = {}
    x, y, z = float(line[31:38]), float(line[39:46]), float(line[47:54])
    element = line[77:78].strip()
    atom['x'], atom['y'], atom['z'], atom['element name'] = x, y, z, element
    return atom


def one_aa_from_pdb(line: str, file: TextIO) \
        -> Tuple[Dict[str, Union[list, int, str]], str]:
    """Loads info about one aminoacid, its atoms
    and coordinates of particular atoms.
    """
    residue_name, residue_sequence_num = line[17:20], line[23:26]
    only_aa = {'residue name': residue_name,
               'residue number': int(residue_sequence_num),
               'atoms': [],
               'coordinates': []}
    while residue_sequence_num == line[23:26]:
        atom_name = line[13:16].strip()
        only_aa['atoms'].append(atom_name)
        only_aa['coordinates'].append(atom_coordinates_pdb(line))
        line = file.readline()
        if line.startswith('TER'):
            # TER -> termination of chain (ignored in this program)
            line = file.readline()
            break
    return only_aa, line


def print_composition_of_protein(aa_counts: Dict[str, int],
                                 most_common: Tuple[str, int]) -> None:
    """Prints count of particular aminoacids in protein.
    Prints the most common aminoacid in protein.
    """
    print(f'Protein consists of:')
    for key, value in aa_counts.items():
        print(f'{value} {key}')
    print()
    print(f'Most common aminoacid is: {most_common[0]} ({most_common[1]} occurrencies).')


def composition_of_protein(aa: List[str]) -> Tuple[Dict[str, int], Tuple[str, int]]:
    """Returns count of particular aminoacids and most common aminoacid."""
    aa_counts = Counter(aa)
    # Counter.most_common returns list of tuples (of x most commons)
    most_common = aa_counts.most_common(1)[0]
    print_composition_of_protein(aa_counts, most_common)
    return aa_counts, most_common


def load_one_standard_aa(file: TextIO) -> Tuple[str, Dict[Tuple[str, str], int]]:
    """Returns name of aminoacid and list of tuples.
    One tuple represents one bond (atom1, atom2, type of bond).
    """
    name = file.readline().strip()
    atoms_of_aa = {}
    line = file.readline()
    # EOF or empty line (end of aa)
    while line and line != '\n':
        first_atom, sec_atom, bond = line.split()
        first_atom, sec_atom = sorted((first_atom, sec_atom))
        atoms_of_aa[(first_atom, sec_atom)] = int(bond)
        line = file.readline()
    return name, atoms_of_aa


def load_standard_aa(path_to_file: str) -> Dict[str, Dict[Tuple[str, str], int]]:
    """Returns standard aminoacids, saved in dictionary.
    {name_of_aa: [(atom1, atom2, type_of_bond), (...), ...]}
    """
    with open(path_to_file) as file:
        aa = {}
        for _ in range(20):
            name, atoms_of_aa = load_one_standard_aa(file)
            aa[name] = atoms_of_aa
        return aa


def read_pdb_file(path_to_file: str) -> None:
    with open(path_to_file) as file:
        line = non_structural_data(file, 'ATOM')
        all_aa_names = []
        all_aa_info = []
        while line.startswith('ATOM'):
            one_aa, line = one_aa_from_pdb(line, file)
            all_aa_names.append(one_aa['residue name'])
            all_aa_info.append(one_aa)
        # print(all_aa_info)
        print()
        composition_of_protein(all_aa_names)
        # load_standard_aa("amino_acids.txt")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='Input file in .sdf or .pdb format')
    args = parser.parse_args()

    input_file = args.file
    if input_file.endswith('.sdf'):
        read_sdf_file(input_file)
    elif input_file.endswith('.pdb'):
        read_pdb_file(input_file)
