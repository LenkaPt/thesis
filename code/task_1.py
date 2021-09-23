'''První úkol.
17.9.2021 - 23.9.2021
Program načítá a zpracovává sdf soubor.
'''

import argparse

def skip_two_lines(file):
    for _ in range(2):
        file.readline()


def fourth_line(file):
    '''Vrací počet atomů a počet vazeb molekuly'''
    skip_two_lines(file)
    line = file.readline()
    return (int(line[:3]), int(line[3:6]))


def print_info_about_molecule(molecule_name, number_of_atoms, number_of_bonds):
    '''Funkce slouží na výpis název sloučeniny, počet atomů vyskytujících se 
    v molekule a počet vazeb.'''
    print(f'Název sloučeniny: {molecule_name}')
    print(f'Počet atomů: {number_of_atoms}.\nPočet vazeb: {number_of_bonds}.\n')


def coordinates(file, number_of_atoms):
    '''Zpracuje část souboru se souřadnicemi atomů.
    Ukládá informace o jednotlivých atomech do slovníku.
    Informace o všech atomech molekuly uloží do seznamu.
    
    Vrací seznam slovníků'''

    coordinates_of_atoms = []
    for _ in range(number_of_atoms):
        atom = {}
        x, y, z, atom_name, *rest = file.readline().split()
        x, y, z = float(x), float(y), float(z)

        atom["x"], atom["y"], atom["z"], atom["atom name"] = x, y, z, atom_name
        coordinates_of_atoms.append(atom)
        
    return coordinates_of_atoms


def bonds(file, number_of_bonds):
    '''Zpracuje část souboru s informacemi o vazbách.
    Uloží informace o každé vazbě do tuplu.
    Informace o vazbách v rámci molekuly jsou uloženy do seznamu.
    
    Vrací seznam tuplů.'''
    
    informations_about_bonds = []
    for i in range(number_of_bonds):
        line = file.readline()
        first_atom, second_atom, type_of_bond = int(line[:3]), int(line[3:6]), int(line[6:9])

        informations_about_bonds.append((first_atom, second_atom, type_of_bond))

    return informations_about_bonds


def molecular_formula(number_of_atoms, molecule_data):
    '''Vypíše sumární vzorec molekuly'''
    atom_names_list = [molecule_data[i]["atom name"] for i in range(number_of_atoms)]
    print("Sumární vzorec: ", end="")
    for atom in (set(atom_names_list)):
        print(f'{atom}{atom_names_list.count(atom)}', end=" ")
    print()


def non_structural_data(file):
    delimiter = "$$$$\n"
    line = file.readline()
    while line != delimiter:
        line = file.readline()


def number_of_bonds_outgoing_of_atom(atoms_bonds_dict, number_of_atom):
    '''Vypíše informaci o počtu vazeb vycházejícího z konkrétního atomu'''
    print(f'Z atomu č. {number_of_atom} vychází {atoms_bonds_dict[number_of_atom][0]}',
    f'jednoduchých vazeb, {atoms_bonds_dict[number_of_atom][1]} dvojných a ',
    f'{atoms_bonds_dict[number_of_atom][2]} trojných vazeb.')


def processed_informations_about_bonds(data_about_bonds, number_of_atoms):
    '''Fce dává dohromady informace o počtu vazeb vycházejících z atomů.
    Vrací slovník, kde klíč je pořadí atomu a hodnota je seznam o třech prvcích.
    Na první pozici je počet jednoduchých vazeb, na druhé počet dvojných vazeb
    a na třetí pořet trojných vazeb vycházejících z daného atomu.'''

    atoms_bonds_dict = {}
    for i in range(1, number_of_atoms +1):
        atoms_bonds_dict[i] = [0, 0, 0]

    for first_atom, sec_atom, type_of_bond in data_about_bonds:
        atoms_bonds_dict[first_atom][type_of_bond - 1] = \
            atoms_bonds_dict.get(first_atom)[type_of_bond - 1] + 1
        atoms_bonds_dict[sec_atom][type_of_bond - 1] = \
            atoms_bonds_dict.get(sec_atom)[type_of_bond - 1] + 1

    number_of_bonds_outgoing_of_atom(atoms_bonds_dict, 2)
    return atoms_bonds_dict


def statistics(number_of_records, names_of_all_molecules_in_file, \
    number_of_atoms_of_each_molecule_in_file):
    '''Vypíše celkový počet molekul v souboru.
    Největší molekulu a nejmenší molekulu.'''
    print(f'Počet záznamů v souboru: {number_of_records}')
    max_atom_count = max(number_of_atoms_of_each_molecule_in_file)
    print(f'Největší počet atomů má molekula: '
        f'{names_of_all_molecules_in_file[number_of_atoms_of_each_molecule_in_file.index(max_atom_count)]}'
        f', konkrétně: {max_atom_count} atomů')
    min_atom_count = min(number_of_atoms_of_each_molecule_in_file)
    print(f'Nejmenší počet atomů má molekula: '
        f'{names_of_all_molecules_in_file[number_of_atoms_of_each_molecule_in_file.index(min_atom_count)]}'
        f', konkrétně: {min_atom_count} atomů')


def read_mol_file(path_to_file):
    with open(path_to_file) as file:
        names_of_all_molecules_in_file = []
        number_of_atoms_of_each_molecule_in_file = []
        number_of_records = 0
        
        while True:
            molecule_name = file.readline().strip()
            if molecule_name == "":
                break

            names_of_all_molecules_in_file.append(molecule_name)
            
            number_of_atoms, number_of_bonds = fourth_line(file)
            number_of_atoms_of_each_molecule_in_file.append(number_of_atoms)

            print_info_about_molecule(molecule_name, number_of_atoms, \
                number_of_bonds)

            coordinates_data = coordinates(file, number_of_atoms)
            data_about_bonds = bonds(file, number_of_bonds)
            processed_informations_about_bonds(data_about_bonds, number_of_atoms)

            molecular_formula(number_of_atoms, coordinates_data)
            non_structural_data(file)

            number_of_records += 1
            print("--------\n")
        statistics(number_of_records, names_of_all_molecules_in_file, \
            number_of_atoms_of_each_molecule_in_file)


parser = argparse.ArgumentParser()
parser.add_argument('file', help="Input file in .sdf format")
args = parser.parse_args()

read_mol_file(args.file)