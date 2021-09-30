"""First task.
17.9.2021 - 23.9.2021
Load and process sdf file.
"""

import argparse
from typing import List, Union, Dict, Tuple, TextIO


def skip_two_lines(file: TextIO) -> None:
    for _ in range(2):
        file.readline()


def fourth_line(file: TextIO) -> (int, int):
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


def coordinates(file: TextIO, number_of_atoms: int) -> List[Dict[str, Union[float, str]]]:
    """Process coordinates part of file.
    Info about one atom is stored in dictionary.
    """
    coordinates_of_atoms = []
    for _ in range(number_of_atoms):
        atom = {}
        x, y, z, atom_name, *rest = file.readline().split()
        x, y, z = float(x), float(y), float(z)

        atom["x"], atom["y"], atom["z"], atom["atom name"] = x, y, z, atom_name
        coordinates_of_atoms.append(atom)

    return coordinates_of_atoms


def bonds(file: TextIO, number_of_bonds: int) -> List[Tuple[int, int, int]]:
    """Process bonds part of file.
    Info about one bond is stored in tuple.
    """
    informations_about_bonds = []
    for i in range(number_of_bonds):
        line = file.readline()
        first_atom, second_atom, type_of_bond = int(line[:3]), int(line[3:6]), int(line[6:9])

        informations_about_bonds.append((first_atom, second_atom, type_of_bond))

    return informations_about_bonds


def molecular_formula(number_of_atoms: int,
                      molecule_data: List[Dict[str, Union[float, str]]]) -> None:
    """Prints molecular formula."""
    atom_names_list = [molecule_data[i]["atom name"] for i in range(number_of_atoms)]
    print("Molecular formula: ", end="")
    for atom in (set(atom_names_list)):
        print(f'{atom}{atom_names_list.count(atom)}', end=" ")
    print()


def non_structural_data(file: TextIO) -> None:
    delimiter = "$$$$\n"
    line = file.readline()
    while line != delimiter:
        line = file.readline()


def bonds_of_atom(atoms_bonds_dict: Dict[int, List[int]],
                  number_of_atom: int) -> None:
    """Prints information about bonds outgoing of specific atom."""
    print(f'Atom number {number_of_atom} forms {atoms_bonds_dict[number_of_atom][0]}',
          f'single bonds, {atoms_bonds_dict[number_of_atom][1]} double bonds and ',
          f'{atoms_bonds_dict[number_of_atom][2]} triple bonds.')


def processed_informations_about_bonds(data_about_bonds: List[Tuple[int, int, int]],
                                       number_of_atoms: int) -> Dict[int, List[int]]:
    """Provides informations about atoms and its bonds.
    Returns dictionary. Key -> number of atom. Value -> list, where
    position means type of bond formed by atom [single, double, triple].
    """
    atoms_bonds_dict = {}
    for i in range(1, number_of_atoms + 1):
        atoms_bonds_dict[i] = [0, 0, 0]

    for first_atom, sec_atom, type_of_bond in data_about_bonds:
        atoms_bonds_dict[first_atom][type_of_bond - 1] = \
            atoms_bonds_dict.get(first_atom)[type_of_bond - 1] + 1
        atoms_bonds_dict[sec_atom][type_of_bond - 1] = \
            atoms_bonds_dict.get(sec_atom)[type_of_bond - 1] + 1

    bonds_of_atom(atoms_bonds_dict, 2)
    return atoms_bonds_dict


def statistics(number_of_records: int, names_of_all_molecules: List[str],
               number_of_atoms_of_each_molecule: List[int]) -> None:
    """Prints out number of all molecules, the largest molecule and
    smallest molecule.
    """
    print(f'Number of molecules: {number_of_records}')
    max_atom_count = max(number_of_atoms_of_each_molecule)
    print(f'Highest number of atoms has molecule: '
          f'{names_of_all_molecules[number_of_atoms_of_each_molecule.index(max_atom_count)]}'
          f', exactly: {max_atom_count} atoms.')
    min_atom_count = min(number_of_atoms_of_each_molecule)
    print(f'Smallest number of atoms has molecule: '
          f'{names_of_all_molecules[number_of_atoms_of_each_molecule.index(min_atom_count)]}'
          f', exactly: {min_atom_count} atoms.')


def read_mol_file(path_to_file: str) -> None:
    with open(path_to_file) as file:
        names_of_all_molecules = []
        number_of_atoms_of_each_molecule = []
        number_of_records = 0

        while True:
            molecule_name = file.readline().strip()
            if molecule_name == "":
                break

            names_of_all_molecules.append(molecule_name)

            number_of_atoms, number_of_bonds = fourth_line(file)
            number_of_atoms_of_each_molecule.append(number_of_atoms)

            print_info_about_molecule(molecule_name, number_of_atoms,
                                      number_of_bonds)

            coordinates_data = coordinates(file, number_of_atoms)
            data_about_bonds = bonds(file, number_of_bonds)
            processed_informations_about_bonds(data_about_bonds, number_of_atoms)

            molecular_formula(number_of_atoms, coordinates_data)
            non_structural_data(file)

            number_of_records += 1
            print("--------\n")
        statistics(number_of_records, names_of_all_molecules,
                   number_of_atoms_of_each_molecule)


parser = argparse.ArgumentParser()
parser.add_argument('file', help="Input file in .sdf format")
args = parser.parse_args()

read_mol_file(args.file)
