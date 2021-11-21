from typing import List, Dict
from collections import Counter


class Atom:
    def __init__(self, x: float, y: float, z: float, name: str, number: int = 0):
        self._x = x
        self._y = y
        self._z = z
        self._name = name
        self._number = number

    @property
    def x(self) -> float:
        return self._x

    @property
    def y(self) -> float:
        return self._y

    @property
    def z(self) -> float:
        return self._z

    @property
    def name(self) -> str:
        return self._name

    @property
    def number(self) -> int:
        return self._number

    def print_coordinates(self) -> None:
        print(f'Coordinates of atom {self._name}: \n\t'
              f'x: {self._x}, y: {self._y}, z: {self._z}')


class Molecule:
    def __init__(self, name: str, atoms: List[Atom], bonds):
        self._name = name
        self._atoms = atoms
        self._bonds = bonds

    @property
    def name(self) -> str:
        return self._name

    @property
    def atoms(self) -> List[Atom]:
        return self._atoms

    def molecular_formula(self) -> str:
        """Returns molecular formula."""
        atoms = [atom.name for atom in self._atoms]
        result = []
        for atom, count in Counter(atoms).items():
            result.append(atom)
            result.append(str(count))

        return ''.join(result)

    def __len__(self) -> int:
        return len(self.atoms)


class Residue:
    def __init__(self, name: str, number: int, atoms: List[Atom]):
        self._name = name
        self._number = number
        self._atoms = atoms

    @property
    def name(self) -> str:
        return self._name

    @property
    def number(self) -> int:
        return self._number

    @property
    def atoms(self) -> List[Atom]:
        return self._atoms

    def get_atom_count(self) -> int:
        return len(self.atoms)

    def __getitem__(self, item: int) -> Atom:
        return self.atoms[item]

    def get_atom_order(self) -> Dict[str, int]:
        """Returns atoms as keys and its order in particular residuum as values"""
        atom_order = {}
        for atom in self:
            atom_order[atom.name] = atom.number

        return atom_order


class Chain:
    def __init__(self, name: str, residues: List[Residue]):
        self._name = name
        self._residues = residues

    @property
    def name(self) -> str:
        return self._name

    @property
    def residues(self) -> List[Residue]:
        return self._residues

    def __len__(self) -> int:
        return len(self.residues)

    def __getitem__(self, item: int) -> Residue:
        return self.residues[item]

    def get_atom_count(self) -> int:
        """Returns number of atoms of particular chain"""
        return sum([residue.get_atom_count() for residue in self])


class Model:
    def __init__(self, name: int, chains: List[Chain]):
        self._name = name
        self._chains = chains

    @property
    def name(self) -> str:
        return self._name

    @property
    def chains(self) -> List[Chain]:
        return self._chains

    def get_number_of_chains(self) -> int:
        return len(self.chains)

    def __getitem__(self, item) -> Chain:
        return self.chains[item]

    def get_residue_count(self) -> int:
        return sum([len(chain) for chain in self])

    def get_most_common_residue(self) -> str:
        residues = []
        for chain in self:
            residues += [residue.name for residue in chain]
        return Counter(residues).most_common(1)[0][0]

    def get_atom_count(self) -> int:
        """Finds out number of atoms in given model of protein."""
        return sum([chain.get_atom_count() for chain in self])


class Protein:
    def __init__(self, name, models: List[Model]):
        self._name = name
        self._models = models

    @property
    def name(self):
        return self._name

    @property
    def models(self) -> List[Model]:
        return self._models

    def get_number_of_models(self) -> int:
        return len(self.models)

    def __getitem__(self, item: int) -> Model:
        return self.models[item]
