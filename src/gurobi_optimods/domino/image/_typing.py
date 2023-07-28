from typing import Any, Protocol, Tuple

import numpy.typing as npt


class Target(Protocol):
    def set_params(self, **kwargs) -> None:
        ...

    def get_targets(self) -> npt.NDArray:
        ...


class ImageHandler(Protocol):
    def __init__(self, img: Any):
        ...

    @property
    def original_size(self) -> Tuple[int, int]:
        ...

    def to_numpy(self, width: int, height: int) -> npt.NDArray:
        ...

    @classmethod
    def from_path(cls, path_to_image: str) -> "ImageHandler":
        ...


class DrawStyle(Protocol):
    @property
    def name(self):
        ...

    @property
    def description(self):
        ...

    @classmethod
    def draw(cls, solution) -> npt.NDArray:
        ...
