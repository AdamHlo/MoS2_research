import numpy as np
import itertools as it
from dataclasses import dataclass
from ase import Atoms, Atom
from ase.io import write
from copy import deepcopy
from typing import Union


@dataclass
class Cell:
    a_1: Union[tuple, list]
    a_2: Union[tuple, list]
    a_3: Union[tuple, list]


@dataclass
class Lattice:
    cell: Cell
    atoms: list

    def get_ase_atoms(self) -> Atoms:
        ase_atoms = Atoms()

        for name, coordinates in self.atoms:
            atom = Atom(symbol=name, position=coordinates)
            ase_atoms.append(atom)

        ase_cell = [list(self.cell.a_1),
                    list(self.cell.a_2),
                    list(self.cell.a_3)]

        ase_atoms.set_pbc(True)
        ase_atoms.set_cell(ase_cell)

        return ase_atoms

    def get_espresso_input(self) -> str:
        lines = 'ATOMIC_POSITIONS angstrom\n'
        for name, coordinates in self.atoms:
            lines += f'{name}    {coordinates[0]}    {coordinates[1]}    {coordinates[2]}\n'

        lines += '\nCELL_PARAMETERS angstrom\n'
        lines += f'{self.cell.a_1[0]}    {self.cell.a_1[1]}    {self.cell.a_1[2]}\n'
        lines += f'{self.cell.a_2[0]}    {self.cell.a_2[1]}    {self.cell.a_2[2]}\n'
        lines += f'{self.cell.a_3[0]}    {self.cell.a_3[1]}    {self.cell.a_3[2]}\n'

        return lines


def generate_lattice(nx, ny, lattice_constant=3.16, sulphur_z_offset=1.58, z_cell_size=50) -> Lattice:
    a_x = np.array([lattice_constant, 0, 0])
    a_y = np.array([-lattice_constant / 2, np.sqrt(3) * lattice_constant / 2, 0])
    centering_vector = np.array([0, 0, z_cell_size / 2])

    s_z_offset = np.array([0, 0, sulphur_z_offset])
    s_y_offset = np.array([0, lattice_constant / np.sqrt(3), 0])

    cell_offset_top = (a_y[1] + s_y_offset[1]) / 2  # arbitrary, just has to be lower than a_y[1] and higher than s_y_offset[1]
    cell_offset_bottom = a_y[1] - cell_offset_top
    pbc_origin_translation = (lattice_constant / 2, cell_offset_bottom, 0)

    atoms = []

    for i, j in it.product(range(nx), range(ny)):
        position = i * a_x + j * a_y + pbc_origin_translation + centering_vector
        atoms.append(('Mo', list(position)))
        atoms.append(('S', list(position + s_y_offset + s_z_offset)))
        atoms.append(('S', list(position + s_y_offset - s_z_offset)))

    e_1 = (nx * lattice_constant, 0, 0)
    e_2 = (-ny * lattice_constant / 2, ny * np.sqrt(3) * lattice_constant / 2, 0)
    e_3 = (0, 0, z_cell_size)

    cell = Cell(e_1, e_2, e_3)

    lattice = Lattice(cell, atoms)

    return lattice


def clone_2d_lattice(lattice: Lattice, n1: int, n2: int) -> Lattice:
    e_1 = np.array(lattice.cell.a_1)
    e_2 = np.array(lattice.cell.a_2)
    e_3 = np.array(lattice.cell.a_3)

    new_atoms = deepcopy(lattice.atoms)

    for i, j in it.product(range(0, n1), range(0, n2)):
        if (i, j) != (0, 0):
            displacement = i * e_1 + j * e_2
            for name, position in lattice.atoms:
                new_position = list(np.array(position) + displacement)
                new_atoms.append((name, new_position))

    e_1_new = n1 * e_1
    e_2_new = n2 * e_2
    e_3_new = e_3

    new_lattice = Lattice(Cell(tuple(e_1_new), tuple(e_2_new), tuple(e_3_new)), new_atoms)

    return new_lattice


def write_ase_configuration(lattices: Union[list[Lattice], Lattice], filename: str = 'output.xyz'):
    if type(lattices) == list:
        atoms_list = []
        for lattice in lattices:
            atoms_list.append(lattice.get_ase_atoms())
        write(filename, atoms_list)
    else:
        atoms = lattices.get_ase_atoms()
        write(filename, atoms)


def main():
    nx = 6
    ny = 6
    clone_n1 = 3
    clone_n2 = 3
    lattice = generate_lattice(nx, ny)
    cloned_lattice = clone_2d_lattice(lattice, clone_n1, clone_n2)
    write_ase_configuration(lattice, 'ase_configuration.xyz')
    write_ase_configuration(cloned_lattice, 'multiplied_ase_configuration.xyz')
    espresso_input = lattice.get_espresso_input()
    print(espresso_input)


if __name__ == '__main__':
    main()
