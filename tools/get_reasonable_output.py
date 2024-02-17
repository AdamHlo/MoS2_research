from xml.dom.minidom import parse, Node
from pathlib import Path
from tools.lattice.create_lattice import Lattice, Cell, clone_2d_lattice, write_ase_configuration
import argparse
import os


BOHR_TO_ANGSTROM = 0.529177249


def get_lattice(file: Path) -> Lattice:
    results = parse(str(file))
    atomic_structure = results.getElementsByTagName('output')[0].getElementsByTagName('atomic_structure')[0]
    atoms_xml = atomic_structure.getElementsByTagName('atom')
    cell_xml = atomic_structure.getElementsByTagName('cell')[0]
    atoms = []

    for at in atoms_xml:
        atom_name = at.attributes.get('name').value
        pos = at.lastChild.nodeValue
        coordinates = [float(value) * BOHR_TO_ANGSTROM for value in pos.strip().split()]
        atoms.append((atom_name, coordinates))

    cell = dict()
    cell_vector_names = ['a1', 'a2', 'a3']

    for name in cell_vector_names:
        cell_vect_str = cell_xml.getElementsByTagName(name)[0].lastChild.nodeValue
        cell_vector = [float(val) * BOHR_TO_ANGSTROM for val in cell_vect_str.strip().split()]
        cell[name] = cell_vector

    lattice_cell = Cell(cell['a1'], cell['a2'], cell['a3'])
    lattice = Lattice(lattice_cell, atoms)

    return lattice


def main():
    print(os.getcwd())
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs='?', default=Path('res/MoS2.xml'), type=Path)
    args = parser.parse_args()
    file = args.file
    lattice = get_lattice(file)
    cloned_lattice = clone_2d_lattice(lattice, 6, 6)
    write_ase_configuration(cloned_lattice, 'multiplied_ase_configuration.xyz')
    espresso_input = lattice.get_espresso_input()
    print(espresso_input)


if __name__ == "__main__":
    main()
