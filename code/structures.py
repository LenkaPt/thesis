from typing import List
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
