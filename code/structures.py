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