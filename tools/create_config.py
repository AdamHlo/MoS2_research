from mako.template import Template
from tools.get_reasonable_output import clone_2d_lattice, get_output_lattice


def generate_config(file='res/mos2.xml'):
    with open('templates/scf.mako') as temp_file:
        template = temp_file.read()

    n1 = 4
    n2 = 4

    output_lattice = get_output_lattice(file)
    # cloned_lattice = clone_2d_lattice(output_lattice, n1, n2)
    atomic_conf = output_lattice.get_espresso_input()
    rendered = Template(template).render(job_name='temp', ecutwfc=40, ecutrho=360, ibrav=0, nat=47, ntyp=2,
                                         nosym='.true.',
                                         conv_thr=0e-6, atomic_configuration=atomic_conf)
    print(rendered)


if __name__ == '__main__':
    generate_config()
