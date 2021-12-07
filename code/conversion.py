import argparse
import sys
from pdb import read_pdb_file
from structures import Protein, Atom
from typing import Dict, Tuple, TextIO
from abc import ABC, abstractmethod
import os


class ProteinToSdf(ABC):
    def __init__(self, protein: Protein, output_file: TextIO):
        self._protein = protein
        self._output_file = output_file

    @property
    def protein(self) -> Protein:
        return self._protein

    @property
    def output_file(self) -> TextIO:
        return self._output_file

    def get_atoms_count_of_protein(self) -> int:
        """Returns number of atoms of protein"""
        # print(f'{self.protein.name}:')
        # i = 1
        # for model in self.protein:
        #     for chain in model:
        #         print(chain.name, chain.get_atom_count())
        #         for residue in chain:
        #             for atom in residue:
        #                 if atom.number != i:
        #                     print(atom.number, atom.name, i)
        #                     return
        #                 i += 1
        # TODO problem with TER
        return self.protein[0].get_atom_count()

    def get_bonds_of_atom(self, atom: Atom, bonds_in_residue: Dict[Tuple[str, str], int],
                          atoms_order_in_residuum: Dict[str, int]) \
            -> Dict[Tuple[str, str], int]:
        """Returns filtered bonds of particular atom within particular residuum."""
        bonds = {}
        for bond, type_of_bond in bonds_in_residue.items():
            first_atom, sec_atom = bond

            if first_atom in atoms_order_in_residuum and sec_atom in atoms_order_in_residuum:
                # bond between two atom is not repeated twice
                if atom.name == first_atom:
                    bonds[bond] = type_of_bond
        return bonds

    def get_bonds_count_of_protein(self) -> int:
        """Returns number of bonds within in the protein"""
        bonds_count = 0
        for chain in self.protein[0]:
            for residue in chain:
                bonds_in_residue = residue.get_standard_bonds()
                atoms_order_in_residue = residue.get_atom_order()
                for atom in residue:
                    bonds_count += len(self.get_bonds_of_atom(atom, bonds_in_residue, atoms_order_in_residue))
                bonds_count += 1  # peptidic bond
        return bonds_count - 1

    def write_name(self) -> None:
        """Writes name of protein into output_file"""
        self.output_file.write(f'{self.protein.name}\n')

    @abstractmethod
    def write_header(self) -> None:
        """Writes sdf header into output_file.
        Particular format implemented in subclass V2000/V3000"""
        pass

    @abstractmethod
    def write_count_line(self) -> None:
        """Writes sdf count_line into output_file.
        Particular format implemented in subclass V2000/V3000"""
        pass

    @abstractmethod
    def write_atom(self, atom: Atom) -> None:
        """Writes atom into sdf output_file.
        Particular format implemented in subclass V2000/V3000"""
        pass

    def write_atom_block(self) -> None:
        """Writes atom block into sdf output_file.
        Particular format implemented in subclass V2000/V3000"""
        # protein[0] == model
        for chain in self.protein[0]:
            for residue in chain:
                for atom in residue:
                    self.write_atom(atom)

    @abstractmethod
    def write_N(self, atoms_order_in_residue) -> None:
        """Writes N of peptidic bond into sdf file.
        Particular format implemented in subclass V2000/V3000"""
        pass

    @abstractmethod
    def write_bonds(self, bonds_of_atom: Dict[Tuple[str, str], int],
                    atoms_order_in_residue: Dict[str, int],
                    bond_number: int) -> None:
        """Writes bonds into sdf output_file.
        Particular format implemented in subclass V2000/V3000"""
        pass

    @abstractmethod
    def write_C(self, atoms_order_in_residue: Dict[str, int],
                bond_number: int) -> None:
        """Writes C of peptidic bond into sdf file.
        Particular format implemented in subclass V2000/V3000"""
        pass

    def write_bond_block(self) -> None:
        """Writes bond block into sdf output file.
        Particular format implemented in subclass V2000/V3000"""
        bond_number = 1
        model = self.protein[0]
        number_of_last_residue = model.get_residue_count()
        for chain in model:
            for residue in chain:
                bonds_in_residue = residue.get_standard_bonds()
                atoms_order_in_residue = residue.get_atom_order()

                if residue.number != 1:
                    self.write_N(atoms_order_in_residue)

                for atom in residue:
                    bonds_of_atom = self.get_bonds_of_atom(atom, bonds_in_residue,
                                                           atoms_order_in_residue)
                    bond_number = self.write_bonds(bonds_of_atom, atoms_order_in_residue, bond_number)

                if residue.number != number_of_last_residue:
                    self.write_C(atoms_order_in_residue, bond_number)
                bond_number += 1    # peptidic bond

    @abstractmethod
    def write_footer(self) -> None:
        """Writes footer of sdf into output file.
        Particular format implemented in subclass V2000/V3000"""
        pass

    @abstractmethod
    def write_to_file(self) -> None:
        """Writes sdf file.
        Particular format implemented in subclass V2000/V3000"""
        pass


class V2000(ProteinToSdf):
    def __init__(self, protein: Protein, output_file: TextIO):
        super().__init__(protein, output_file)

    def write_header(self):
        self.output_file.write('\n\n')

    def write_count_line(self, atom_lists: int = 0, fff: int = 0,
                         chiral_flag: int = 0, stext_entries: int = 0,
                         xxx: int = 0, rrr: int = 0, ppp: int = 0,
                         iii: int = 0, properties_line: int = 999,
                         v: str = 'V2000'):
        self.output_file.write(f'{self.get_atoms_count_of_protein():>3}'
                               f'{self.get_bonds_count_of_protein():>3}'
                               f'{atom_lists:>3}'
                               f'{fff:>3}'
                               f'{chiral_flag:>3}'
                               f'{stext_entries:>3}'
                               f'{xxx:>3}'
                               f'{rrr:>3}'
                               f'{ppp:>3}'
                               f'{iii:>3}'
                               f'{properties_line:>3}'
                               f'{v:>6}'
                               f'\n')

    def write_atom(self, atom: Atom, mass_difference: int = 0, charge: int = 0,
                   atom_stereo_parity: int = 0, hydrogen_count: int = 0,
                   stereo: int = 0, valence: int = 0, ho_designator: int = 0,
                   rrr: int = 0, iii: int = 0, mapping_number: int = 0,
                   inversion_fag: int = 0, change_flag: int = 0):
        self.output_file.write(f'{atom.x:>10.4f}'
                               f'{atom.y:>10.4f}'
                               f'{atom.z:>10.4f}'
                               f' '
                               f'{atom.name[0]:<3}'
                               f'{mass_difference:>2}'
                               f'{charge:>3}'
                               f'{atom_stereo_parity:>3}'
                               f'{hydrogen_count:>3}'
                               f'{stereo:>3}'
                               f'{valence:>3}'
                               f'{ho_designator:>3}'
                               f'{rrr:>3}{iii:>3}'
                               f'{mapping_number:>3}'
                               f'{inversion_fag:>3}'
                               f'{change_flag:>3}'
                               f'\n')

    def write_N(self, atoms_order_in_residue: Dict[str, int]):
        self.output_file.write(f'{atoms_order_in_residue["N"]:>3}'
                               f'{1:>3}'
                               f'{0:>3}'
                               f'{0:>3}'
                               f'{0:>3}'
                               f'{0:>3}'
                               f'\n')

    def write_bonds(self, bonds_of_atom: Dict[Tuple[str, str], int],
                    atoms_order_in_residue: Dict[str, int],
                    bond_number: int,
                    bond_stereo: int = 0, xxx: int = 0,
                    bond_topology: int = 0,
                    reacting_center_status: int = 0) -> int:
        """Function writes bonds of particular atom into output file (sdf V2000 format)"""
        for bond, type_of_bond in bonds_of_atom.items():
            first_atom, sec_atom = bond

            self.output_file.write(f'{atoms_order_in_residue[first_atom]:>3}'
                                   f'{atoms_order_in_residue[sec_atom]:>3}'
                                   f'{type_of_bond:>3}'
                                   f'{bond_stereo:>3}'
                                   f'{xxx:>3}'
                                   f'{bond_topology:>3}'
                                   f'{reacting_center_status:>3}'
                                   f'\n')
            bond_number += 1
        return bond_number

    def write_C(self, atoms_order_in_residue: Dict[str, int], bond_number: int):
        self.output_file.write(f'{atoms_order_in_residue["C"]:>3}')

    def write_footer(self):
        self.output_file.write(f'M  END'
                               f'\n'
                               f'$$$$'
                               f'\n')

    def write_to_file(self) -> None:
        self.write_name()
        self.write_header()
        self.write_count_line()
        self.write_atom_block()
        self.write_bond_block()
        self.write_footer()


class V3000(ProteinToSdf):
    def __init__(self, protein: Protein, output_file: TextIO):
        super().__init__(protein, output_file)

    def write_header(self):
        self.output_file.write(f'\n\n'
                               f' 0 0'
                               f'{999:>29}'
                               f'{"V3000":>6}'
                               f'\n')

    def write_count_line(self, sgroups: int = 0,
                         num_of_3d_constraints: int = 0,
                         chiral: int = 0):
        self.output_file.write(f'M  V30 COUNTS '
                               f'{self.get_atoms_count_of_protein()} '
                               f'{self.get_bonds_count_of_protein()} '
                               f'{sgroups} '
                               f'{num_of_3d_constraints} '
                               f'{chiral}'
                               f'\n')

    def write_atom(self, atom: Atom):
        self.output_file.write(f'M  V30 '
                               f'{atom.number} '
                               f'{atom.name[0]} '
                               f'{atom.x} '
                               f'{atom.y} '
                               f'{atom.z}'
                               f'\n')

    def write_N(self, atoms_order_in_residue: Dict[str, int]):
        self.output_file.write(f'{atoms_order_in_residue["N"]}\n')

    def write_bonds(self, bonds_of_atom: Dict[Tuple[str, str], int],
                    atoms_order_in_residue: Dict[str, int],
                    bond_number: int) -> int:
        for bond, type_of_bond in bonds_of_atom.items():
            first_atom, sec_atom = bond

            self.output_file.write(f'M  V30 {bond_number} {type_of_bond} '
                                   f'{atoms_order_in_residue[first_atom]} '
                                   f'{atoms_order_in_residue[sec_atom]}\n')
            bond_number += 1
        return bond_number

    def write_C(self, atoms_order_in_residue: Dict[str, int], bond_number: int):
        self.output_file.write(f'M  V30 {bond_number} {1} {atoms_order_in_residue["C"]} ')

    def write_footer(self):
        self.output_file.write(f'M  END\n')

    def write_to_file(self):
        self.write_name()
        self.write_header()

        self.output_file.write(f'M  V30 BEGIN CTAB\n')
        self.write_count_line()

        self.output_file.write(f'M  V30 BEGIN ATOM\n')
        self.write_atom_block()
        self.output_file.write(f'M  V30 END ATOM\n')

        self.output_file.write(f'M  V30 BEGIN BOND\n')
        self.write_bond_block()
        self.output_file.write(f'M  V30 END BOND\n')
        self.output_file.write(f'M  V30 END CTAB\n')

        self.write_footer()


def convert_pdb_to_sdf(protein: Protein, destination_file):
    with open(destination_file, mode='a', encoding='utf8') as output_file:
        if protein[0].get_atom_count() <= 999:
            data = V2000(protein, output_file)
        else:
            data = V3000(protein, output_file)
        data.write_to_file()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='Input file in .sdf or .pdb format')
    parser.add_argument('-c', '--conversion', dest='conversion_file',
                        help='Allows to convert protein from input pdb file to output sdf file. '
                             'Please specify name of output .sdf file.',
                        default=False)
    args = parser.parse_args()

    input_file = args.file

    try:
        protein = read_pdb_file(input_file)
    except ValueError as e:
        print(e)
        sys.exit(1)

    if args.conversion_file:
        convert_pdb_to_sdf(protein, args.conversion_file)

        os.startfile(args.conversion_file)
