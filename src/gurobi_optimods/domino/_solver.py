"""
Classes:
---------
- Solver
"""
import itertools
from typing import Callable, Dict, List, Optional, Tuple

import gurobipy as gp
import numpy as np
import numpy.typing as npt

from gurobi_optimods import _data_dir
from gurobi_optimods.domino._structs import Domino, Placement, Position, Solution
from gurobi_optimods.domino.image._base import (
    TargetArray,
    TargetImage,
    get_image_handler,
    matplotlib_array_plot,
)
from gurobi_optimods.domino.image._typing import Target
from gurobi_optimods.utils import optimod

_example_image_filepath = str(_data_dir.parent / "data" / "domino" / "keanu.jpg")


def _dominoes_in_set(max_dots):
    """How many dominoes in a double-x set?"""
    return (max_dots + 1) * (max_dots + 2) / 2


class Problem:
    """A class used for encapsulating problem data.  Used for caching on Solver class.
    Not exposed to user."""

    def __init__(self, targets: npt.NDArray):
        self.height, self.width = targets.shape
        self.max_dots: int = int(targets.max())
        self.positions: List[Position] = self._create_positions()
        self.dominoes: List[Domino] = self._create_dominoes()
        self.placements: Dict[
            Tuple[Domino, Position], Placement
        ] = self._create_placements(targets)
        self.grid_to_position_map: Dict[
            int, List[Position]
        ] = self._map_grid_to_position()

    def _create_positions(self) -> List[Position]:
        # Will need modification for non-rectangle arrangements
        positions = [
            Position(r, c, True)
            for r in range(self.height)
            for c in range(self.width - 1)
        ]
        positions.extend(
            [
                Position(r, c, False)
                for r in range(self.height - 1)
                for c in range(self.width)
            ]
        )
        return positions

    def _create_dominoes(self) -> List[Domino]:
        return list(
            Domino(n1, n2)
            for n1, n2 in itertools.combinations_with_replacement(
                range(self.max_dots + 1), 2
            )
        )

    def _create_placements(
        self, targets: npt.NDArray
    ) -> Dict[Tuple[Domino, Position], Placement]:
        def create_placement(domino, position):
            placement = Placement(domino, position)
            self._calculate_placement_objective(targets, placement)
            return placement

        return {
            (domino, position): create_placement(domino, position)
            for domino in self.dominoes
            for position in self.positions
        }

    def _calculate_placement_objective(self, targets, placement: Placement) -> None:
        """Given a Placement object, calculates the best orientation of the domino
        and sets this orientation and best objective on the Placement object.
        """

        def calc_objective(target: npt.NDArray, domino: Domino):
            # this works with numpy.sub since Domino is a named tuple
            return np.square(target - domino).sum()

        target = self._get_array_values(targets, placement.position)
        objective = calc_objective(target, placement.domino)
        flipped_domino = placement.domino.flip()
        flipped_objective = calc_objective(target, flipped_domino)
        if objective <= flipped_objective:
            placement.objective = objective
            placement.placed_domino = placement.domino
        else:
            placement.objective = flipped_objective
            placement.placed_domino = flipped_domino

    def _get_array_values(self, array: npt.NDArray, position: Position) -> np.ndarray:
        """generic helper function which return a 1D numpy array of length 2
        taken from the values in the array argument"""
        if position.horizontal:
            return array[position.row, position.col : position.col + 2]
        return array[position.row : position.row + 2, position.col]

    def _map_grid_to_position(self) -> Dict[int, List[Position]]:
        """Used to ensure dominoes overlap and no squares are uncovered"""
        grid_indices = np.array(range(self.height * self.width)).reshape(
            (self.height, self.width)
        )

        grid_to_position_map: Dict[int, List[Position]] = {}
        for position in self.positions:
            for grid_index in tuple(self._get_array_values(grid_indices, position)):
                if grid_index not in grid_to_position_map:
                    grid_to_position_map[grid_index] = []
                grid_to_position_map[grid_index].append(position)

        return grid_to_position_map

    @property
    def min_domino_sets_required(self) -> int:
        dominoes_in_set = _dominoes_in_set(self.max_dots)
        return int(np.ceil(self.height * self.width / 2 / dominoes_in_set))


class Solver:
    """
    A class used to solve the Domino Art problem.
    The only class which is intended to be exposed to user.
    """

    def __init__(
        self, target_obj: Target, array_plotter: Callable = matplotlib_array_plot
    ):
        self._target_obj: Target = target_obj
        self._array_plotter: Callable = array_plotter
        self._problem: Optional[Problem] = None

    def set_params(self, **kwargs) -> None:
        """Sets parameters which are required for translating an image into
        an array of target values which will be matched as closely as possible
        by dominoes.  Will show target image if show=True.

        Parameters
        ------
        **kwargs:
            Keyword parameters for
            - width (int)
            - max_dots (int)
            - show (bool), default True
        """
        show = kwargs.pop("show", True)
        self._target_obj.set_params(**kwargs)
        self._problem = None  # reset cache
        if show:
            self.plot_targets()

    def plot_targets(self) -> None:
        """Renders the target array with grayscale values"""
        self._array_plotter(self._target_obj.get_targets())

    @optimod()
    def find_optimal_placements(
        self, *, max_domino_sets: int, create_env: gp.Env
    ) -> Solution:
        """Builds and solves optimization model to place the dominoes
        in a configuration to most accurately approximate the target image.

        Parameters
        ----------
        max_domino_sets : int
            The number of sets of dominos to use.  Determines the maximum amount of
            times a domino tile can be replicated.

        Returns
        -------
        Solution

        Raises
        -------
        ValueError
            Error raised when the the number of dominoes (a function of max_dots
            and max_domino_sets) is insufficient to cover the target array.
        """
        if self._problem is None:  # cache check
            targets = self._target_obj.get_targets()
            self._problem = Problem(targets)

        if max_domino_sets < self._problem.min_domino_sets_required:
            raise ValueError(
                f"""The argument for max_domino_sets needs to be at least
                {self._problem.min_domino_sets_required} given the size of the image"""
            )

        with create_env() as env, gp.Model(env=env) as model:
            model = gp.Model()
            x = model.addVars(
                self._problem.dominoes, self._problem.positions, vtype="B"
            )

            # optimal arrangement and objective of placements already calculated
            model.setObjective(
                gp.quicksum(
                    x[d, p] * self._problem.placements[d, p].objective
                    for p in self._problem.positions
                    for d in self._problem.dominoes
                )
            )

            # number of times each domino is replicated must be no greater than
            #  the maximum number of domino sets used
            model.addConstrs(
                gp.quicksum(x[d, p] for p in self._problem.positions) <= max_domino_sets
                for d in self._problem.dominoes
            )

            # no overlaps, no uncovered squares
            for positions in self._problem.grid_to_position_map.values():
                model.addConstr(
                    gp.quicksum(
                        x[d, p] for d in self._problem.dominoes for p in positions
                    )
                    == 1
                )

            model.optimize()

            return Solution(
                width=self._problem.width,
                height=self._problem.height,
                max_dots=self._problem.max_dots,
                num_dominoes_in_set=_dominoes_in_set(self._problem.max_dots),
                placements=tuple(
                    placement
                    for placement in self._problem.placements.values()
                    if x[placement.domino, placement.position].X >= 0.99
                ),
            )

    @classmethod
    def from_target_array(cls, *, array: npt.NDArray) -> "Solver":
        if not isinstance(array, np.ndarray):
            raise ValueError(
                f"Target array is of type {type(array)} but should be a numpy array."
            )
        if not array.ndim == 2:
            raise ValueError(
                f"Target array has {array.ndim} dimensions but should have only 2."
            )
        if array.shape[0] % 2 and array.shape[1] % 2:
            raise ValueError(
                "Target array requires resizing.  The value of width*height must divisible by 2."
            )
        if array.min() < 0:
            raise ValueError(
                f"Target array contains value {array.min()} but should only have values between 0 and 9 (inclusive)."
            )
        if array.max() > 9:
            raise ValueError(
                f"Target array contains value {array.max()} but should only have values between 0 and 9 (inclusive)."
            )
        if (np.round(array) - array).any():
            raise ValueError("Target array must only contain integer values.")
        return cls(TargetArray(array))

    @classmethod
    def from_file(cls, *, filepath: str, backend="pillow") -> "Solver":
        """Creates a new Solver instance from an "image" object

        Parameters
        ----------
        filepath: str
            Location of an image file, such as a jpeg or png image.
        backend : {"pillow", "opencv", "skimage"}
            Indicates the choice of image handler for resizing, grayscale conversions
            and conversion from an image format to numpy array.

        Returns
        -------
        Solver
        """
        image_handler_cls = get_image_handler(backend)
        return cls(TargetImage(image_handler_cls.from_path(filepath)))

    @classmethod
    def from_image(cls, *, image, backend="pillow") -> "Solver":
        """Creates a new Solver instance from an "image" object

        Parameters
        ----------
        image:
            This is an image object, as determined by the image handler backend.
            For "pillow" this is an instance of PIL.Image.
            For "opencv" and "skimage" this is a numpy array.
        backend : {"pillow", "opencv", "skimage"}
            Indicates the choice of image handler for resizing, grayscale conversions
            and conversion from an image format to numpy array.

        Returns
        -------
        Solver
        """
        image_handler_cls = get_image_handler(backend)
        return cls(TargetImage(image_handler_cls(image)))

    @classmethod
    def from_example(cls, *, backend="pillow") -> "Solver":
        """Creates a new Solver instance from an example image

        Parameters
        ----------
        backend : {"pillow", "opencv", "skimage"}
            Indicates the choice of image handler for resizing, grayscale conversions
            and conversion from an image format to numpy array.

        Returns
        -------
        Solver
        """
        return cls.from_file(filepath=_example_image_filepath, backend=backend)
