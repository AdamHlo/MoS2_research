from xml.dom.minidom import parse, Node


BOHR_TO_ANGSTROM = 0.529177249


if __name__ == "__main__":
    path = 'res/MoS2.xml'
    results = parse(path)
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

