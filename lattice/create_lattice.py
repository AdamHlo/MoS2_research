import numpy as np
import itertools as it
from dataclasses import dataclass
from ase import Atoms, Atom
from ase.io import write


@dataclass
class Cell:
    e_1: tuple
    e_2: tuple
    e_3: tuple


@dataclass
class Lattice:
    cell: Cell
    Mo_config: np.array
    S_config_up: np.array
    S_config_down: np.array

    def get_ase_atoms(self) -> Atoms:
        atoms = Atoms()

        def add_ase_atoms(symbol, configuration):
            for i in range(configuration.shape[0]):
                coordinates = configuration[i, :]
                atom = Atom(symbol=symbol, position=tuple(coordinates))
                atoms.append(atom)

        add_ase_atoms('Mo', self.Mo_config)
        add_ase_atoms('S', self.S_config_down)
        add_ase_atoms('S', self.S_config_up)

        ase_cell = [list(self.cell.e_1),
                    list(self.cell.e_2),
                    list(self.cell.e_3)]

        atoms.set_pbc(True)
        atoms.set_cell(ase_cell)

        return atoms


def generate_lattice(nx, ny, lattice_constant=3.16, sulphur_z_offset=1.58, z_cell_size=50) -> Lattice:
    a_x = np.array([lattice_constant, 0, 0])
    a_y = np.array([-lattice_constant / 2, np.sqrt(3) * lattice_constant / 2, 0])
    centering_vector = np.array([0, 0, z_cell_size / 2])

    s_z_offset = np.array([0, 0, sulphur_z_offset])
    s_y_offset = np.array([0, lattice_constant / np.sqrt(3), 0])

    cell_offset_top = (a_y[1] + s_y_offset[1]) / 2  # arbitrary, just has to be lower than a_y[1] and higher than s_y_offset[1]
    cell_offset_bottom = a_y[1] - cell_offset_top
    pbc_origin_translation = (lattice_constant / 2, cell_offset_bottom, 0)

    Mo_config = []
    S_config_up = []
    S_config_down = []

    for i, j in it.product(range(nx), range(ny)):
        position = i * a_x + j * a_y + pbc_origin_translation + centering_vector
        Mo_config.append(position)
        S_config_up.append(position + s_y_offset + s_z_offset)
        S_config_down.append(position + s_y_offset - s_z_offset)

    e_1 = (nx * lattice_constant, 0, 0)
    e_2 = (-ny * lattice_constant / 2, ny * np.sqrt(3) * lattice_constant / 2, 0)
    e_3 = (0, 0, z_cell_size)

    cell = Cell(e_1, e_2, e_3)

    lattice = Lattice(cell, np.stack(Mo_config, axis=0), np.stack(S_config_up, axis=0), np.stack(S_config_down, axis=0))

    return lattice


def clone_2d_lattice(lattice: Lattice, order: int = 1) -> Lattice:
    Mo_config = lattice.Mo_config
    S_config_up = lattice.S_config_up
    S_config_down = lattice.S_config_down

    e_1 = np.array(lattice.cell.e_1)
    e_2 = np.array(lattice.cell.e_2)
    e_3 = np.array(lattice.cell.e_3)

    origin_displacement = order * (e_1 + e_2)

    for i, j in it.product(range(-order, order+1), range(-order, order+1)):
        if (i, j) != (0, 0):
            displacement = i * e_1 + j * e_2
            Mo_config = np.append(Mo_config, lattice.Mo_config + displacement, axis=0)
            S_config_up = np.append(S_config_up, lattice.S_config_up + displacement, axis=0)
            S_config_down = np.append(S_config_down, lattice.S_config_down + displacement, axis=0)

    Mo_config = Mo_config + origin_displacement
    S_config_up = S_config_up + origin_displacement
    S_config_down = S_config_down + origin_displacement

    e_1_new = (2 * order + 1) * e_1
    e_2_new = (2 * order + 1) * e_2
    e_3_new = e_3

    new_lattice = Lattice(Cell(tuple(e_1_new), tuple(e_2_new), tuple(e_3_new)), Mo_config, S_config_up, S_config_down)

    return new_lattice


def write_ase_configuration(nx: int, ny: int, filename: str = 'output.xyz'):
    lattice = generate_lattice(nx, ny)
    atoms = lattice.get_ase_atoms()
    write(filename, atoms)


def write_ase_configuration_extended(nx: int, ny: int, filename: str = 'output.xyz'):
    """
    Used to visually check whether cell for PBCs is properly set (alignment of multiplied lattice)
    :param nx:
    :param ny:
    :param filename:
    :return:
    """
    lattice = generate_lattice(nx, ny)
    extended_lattice = clone_2d_lattice(lattice)
    atoms = extended_lattice.get_ase_atoms()
    write(filename, atoms)


def main():
    write_ase_configuration(6, 6, 'ase_configuration.xyz')
    write_ase_configuration_extended(6, 6, 'multiplied_ase_configuration.xyz')


if __name__ == '__main__':
    main()
