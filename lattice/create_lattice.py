import numpy as np
import itertools as it
from dataclasses import dataclass


@dataclass
class CellXY:
    top_left: tuple
    bottom_left: tuple
    top_right: tuple
    bottom_right: tuple


@dataclass
class Cell:
    cell_xy: CellXY
    z_bottom: float
    z_top: float


@dataclass
class Result:
    cell: Cell
    Mo_config: list[np.array]
    S_config_up: list[np.array]
    S_config_down: list[np.array]


def generate_triangular_Mo(lattice_constant, nx, ny):
    a_x = np.array([lattice_constant, 0, 0])
    a_y = np.array([lattice_constant / 2, np.sqrt(3) * lattice_constant / 2, 0])

    Mo_config = []

    for i, j in it.product(range(nx), range(ny)):
        position = i * a_x + j * a_y
        Mo_config.append(position)

    return Mo_config
