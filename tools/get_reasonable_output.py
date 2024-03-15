from xml.dom.minidom import parse, Node
from pathlib import Path
from tools.lattice.create_lattice import Lattice, Cell, clone_2d_lattice, write_ase_configuration
import argparse


BOHR_TO_ANGSTROM = 0.529177249


def get_output_lattice(file: Path) -> Lattice:
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


def get_sequence_of_lattice_steps(file: Path) -> list[Lattice]:
    results = parse(str(file))
    steps = results.getElementsByTagName('step')

    lattices = []

    for step in steps:
        atoms = []
        atomic_structure = step.getElementsByTagName('atomic_structure')[0]
        atoms_xml = atomic_structure.getElementsByTagName('atom')
        cell_xml = atomic_structure.getElementsByTagName('cell')[0]

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
        lattices.append(lattice)

    return lattices


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs='?', default=Path('res/mos2_bomd_temp_rescale.xml'), type=Path)
    args = parser.parse_args()
    file = args.file
    n1 = 4
    n2 = 4

    lattice_steps = get_sequence_of_lattice_steps(file)
    output_lattice = get_output_lattice(file)
    lattice_steps.append(output_lattice)

    cloned_lattice = clone_2d_lattice(output_lattice, n1, n2)
    write_ase_configuration(output_lattice, 'output_configuration.xyz')
    write_ase_configuration(lattice_steps, 'sequence.xyz')
    write_ase_configuration(cloned_lattice, 'cloned_configuration.xyz')

    cloned_lattice_espresso = cloned_lattice.get_espresso_input()
    output_lattice_espresso = output_lattice.get_espresso_input()
    print(f"----------------- cloned lattice |{n1}x{n2}| -------------------\n\n")
    print(cloned_lattice_espresso)
    print("\n\n")
    print(f"----------------------- output lattice -------------------------\n\n")
    print(output_lattice_espresso)


if __name__ == "__main__":
    main()
