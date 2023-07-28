from collections import Counter
from dataclasses import dataclass
from typing import NamedTuple, Optional, Tuple


class Domino(NamedTuple):
    num1: int
    num2: int

    def flip(self):
        return type(self)(self.num2, self.num1)


class Position(NamedTuple):
    row: int
    col: int
    horizontal: bool


@dataclass
class Placement:
    domino: Domino
    position: Position
    _placed_domino: Optional[Domino] = None
    objective: Optional[float] = None

    @property
    def placed_domino(self) -> Domino:
        if self._placed_domino is None:
            raise RuntimeError("Problem has not been setup correctly")
        return self._placed_domino

    @placed_domino.setter
    def placed_domino(self, domino: Domino) -> None:
        self._placed_domino = domino


@dataclass
class Solution:
    """
    A class used to represent a solution to the Domino Art problem
    """

    width: int
    height: int
    max_dots: int
    num_dominoes_in_set: int
    placements: Tuple[Placement, ...]

    def print_stats(self):
        """Prints to screen several solution statistics"""
        num_dominoes_used = int(self.width * self.height / 2)
        num_sets_used = max(
            Counter([placement.domino for placement in self.placements]).values()
        )
        total_dominoes_in_sets = num_sets_used * self.num_dominoes_in_set

        print(
            f"""
Solution Stats
----------------------------
Width in squares:  {self.width}
Height in squares: {self.height}
Dominoes used:     {num_dominoes_used}
Domino sets used:  {num_sets_used}
Dominoes used (%): {num_dominoes_used/total_dominoes_in_sets:.1%}
        """
        )
