from mako.template import Template
from tools.get_reasonable_output import clone_2d_lattice, get_output_lattice
from tools.lattice.create_lattice import Lattice
from dataclasses import dataclass, asdict
from typing import Union
from pathlib import Path


def bool2espresso(value: bool):
    if value:
        return '.true.'
    else:
        return '.false.'


def float2espresso(value: float):
    return '{:f}'.format(value)


@dataclass
class ScfParameters:
    job_name: str
    ecutwfc: int = 50
    ecutrho: int = 400
    ibrav: int = 0
    nosym: str = bool2espresso(True)
    conv_thr: str = float2espresso(1e-6)  # TODO: potential source of mistakes, handle it systematically in the future


@dataclass
class CPParameters:
    job_name: str
    ecutwfc: int = 50
    ecutrho: int = 400
    ibrav: int = 0
    conv_thr: str = float2espresso(1e-6)  # TODO: same as above
    ndr: int = 50
    ndw: int = 50


@dataclass
class Parameters:
    source_lattice_file: Path
    template_file: Path
    espresso_output_file: Path
    simulation_parameters: Union[ScfParameters, CPParameters]


def generate_espresso_config(simulation_parameters,
                             nat: int,
                             ntyp: int,
                             lattice: Lattice,
                             template_file: Path) -> str:
    atomic_conf = lattice.get_espresso_input()

    with open(template_file) as file:
        template_text = file.read()

    template = Template(template_text)
    rendered = template.render(**asdict(simulation_parameters),
                               ntyp=ntyp,
                               nat=nat,
                               atomic_configuration=atomic_conf)
    return rendered


def render_config(parameters: Parameters, output_lattice: Lattice):

    nat = output_lattice.size
    ntyp = 2
    rendered_config = generate_espresso_config(parameters.simulation_parameters, nat, ntyp, output_lattice,
                                               parameters.template_file)
    file_path = (parameters.espresso_output_file / parameters.simulation_parameters.job_name).with_suffix('.sh')
    with open(file_path, mode='w+') as file:
        file.write(rendered_config)


def generate_cp_pw_test_configs():
    for i in range(4):
        randomize_scale = 0.05

        pw_parameters = Parameters(
            source_lattice_file=Path('res/mos2.xml'),
            template_file=Path('templates/pw_scf.mako'),
            espresso_output_file=Path('../pw/'),
            simulation_parameters=ScfParameters(
                job_name=f'4x4_s_vacancy_randomized_lattice_scf_{i}'
            )
        )

        cp_parameters = Parameters(
            source_lattice_file=Path('res/mos2.xml'),
            template_file=Path('templates/cp_scf.mako'),
            espresso_output_file=Path('../cp/'),
            simulation_parameters=CPParameters(
                job_name=f'4x4_s_vacancy_randomized_lattice_cp_scf_{i}'
            )
        )

        output_lattice = get_output_lattice(pw_parameters.source_lattice_file).noised_lattice(randomize_scale)

        render_config(pw_parameters, output_lattice)
        render_config(cp_parameters, output_lattice)


def generate_cp_dynamics_config():
    pass


def main():
    generate_cp_pw_test_configs()


if __name__ == '__main__':
    main()
