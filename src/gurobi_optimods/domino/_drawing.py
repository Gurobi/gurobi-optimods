import itertools
import warnings
from typing import Dict, List

import numpy as np
import numpy.typing as npt

from gurobi_optimods.domino._solver import Solution
from gurobi_optimods.domino._structs import Placement
from gurobi_optimods.domino.image._typing import DrawStyle
from gurobi_optimods.utils import run


class BasicDotStyle:
    """implements protocol DrawStyle"""

    black_dots: bool

    @run
    @staticmethod
    def domino_squares():
        dot = np.ones((5, 5))
        dot[[0, 0, 4, 4], [0, 4, 0, 4]] = 0

        # upper left corners
        dot_origins = dict(enumerate(itertools.product((1, 8, 15), (1, 8, 15))))

        def make_domino_square(*args):
            """Given indices of dots, creates a 'drawing' of the domino square
            as a binary numpy array.

            Index of dots on a 9-dot tile:
            +-------+
            | 0 1 2 |
            | 3 4 5 |
            | 6 7 8 |
            +-------+

            Returns
            -------
            numpy.ndarray
            """
            square = np.zeros((21, 21))
            for dot_key in args:
                origin_y, origin_x = dot_origins[dot_key]
                square[origin_y : origin_y + 5, origin_x : origin_x + 5] = dot
            return square

        return {
            0: make_domino_square(),
            1: make_domino_square(4),
            2: make_domino_square(2, 6),
            3: make_domino_square(2, 4, 6),
            4: make_domino_square(0, 2, 6, 8),
            5: make_domino_square(0, 2, 4, 6, 8),
            6: make_domino_square(0, 1, 2, 6, 7, 8),
            7: make_domino_square(0, 1, 2, 4, 6, 7, 8),
            8: make_domino_square(0, 1, 2, 3, 5, 6, 7, 8),
            9: make_domino_square(0, 1, 2, 3, 4, 5, 6, 7, 8),
        }

    @classmethod
    def draw(cls, solution: Solution) -> npt.NDArray:
        """Creates an image made from placing dominoes according to a Solution.

        Parameters
        ----------
        solution : Solution
            As returned by Problem.find_optimal_placements

        Returns
        ----------
        numpy.ndarray
        """
        return cls._draw(solution)

    @classmethod
    def _draw(cls, solution: Solution) -> npt.NDArray:
        def make_tile(placement: Placement):
            tile = np.zeros((23, 23 * 2))
            num1 = placement.placed_domino.num1
            num2 = placement.placed_domino.num2

            if not cls.black_dots:
                num1 = solution.max_dots - num1
                num2 = solution.max_dots - num2

            tile[1:22, 1:22] = cls.domino_squares[num1] * 2  # first square
            tile[1:22, 24:45] = cls.domino_squares[num2] * 2  # second square

            # domino border
            tile[0, :] = 1
            tile[-1, :] = 1
            tile[:, 0] = 1
            tile[:, -1] = 1

            if not placement.position.horizontal:
                tile = np.rot90(tile, -1)
            return tile

        canvas = np.ones((23 * solution.height, 23 * solution.width))
        for placement in solution.placements:
            tile = make_tile(placement)
            tile_px_height, tile_px_width = tile.shape
            origin_y, origin_x = (
                placement.position.row * 23,
                placement.position.col * 23,
            )
            canvas[
                origin_y : origin_y + tile_px_height,
                origin_x : origin_x + tile_px_width,
            ] = tile

        if not cls.black_dots:
            canvas = canvas.max() - canvas

        return canvas


class WhiteTileDotStyle(BasicDotStyle):
    black_dots = True

    @property
    def name(self):
        return "black_dots"

    @property
    def description(self):
        return "black dots on white tiles"


class BlackTileDotStyle(BasicDotStyle):
    black_dots = False

    @property
    def name(self):
        return "white_dots"

    @property
    def description(self):
        return "white dots on black tiles"


styles: List[DrawStyle] = [
    WhiteTileDotStyle(),
    BlackTileDotStyle(),
]

style_map: Dict[str, DrawStyle] = {style.name: style for style in styles}
style_keys = set(style_map.keys())


def draw(*, style, solution: Solution) -> npt.NDArray:
    """
    Draws the solution according to a style

    Parameters
    ----------
    style : str
        Choose from the following styles
            * {style_descriptions}
    solution : Solution
        A solution object as created by the Solver class

    Returns
    -------
    np.ndarray

    Raises
    ------
    ValueError
        Raised whenever style_name is not a valid choice
    """
    if style not in style_keys:
        raise ValueError(f"Argument for style_name must be one of {style_keys}")
    array = style_map[style].draw(solution)
    try:
        import matplotlib.pyplot as plt

        ax = plt.imshow(array, cmap="Greys")
        ax.axes.set_axis_off()
    except:
        warnings.warn("Matplotlib must be installed to draw the solution")
    return array


def style_descriptions(sep):
    return f"{sep}".join(
        [" - ".join([f'"{style.name}"', style.description]) for style in styles]
    )


assert draw.__doc__ is not None
draw.__doc__ = draw.__doc__.format(
    style_descriptions=style_descriptions("\n                * ")
)
