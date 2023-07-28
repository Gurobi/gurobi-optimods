import warnings
from typing import Optional

import numpy as np
import numpy.typing as npt

from gurobi_optimods.domino.image._typing import ImageHandler

try:
    import matplotlib.pyplot as plt

    has_matplotlib = True
except ImportError:
    has_matplotlib = False


def matplotlib_array_plot(array: npt.NDArray) -> None:
    if not has_matplotlib:
        raise ImportError(f"Matplotlib is required to show images.")
    ax = plt.imshow(array, cmap="Greys")
    ax.axes.set_axis_off()


def get_image_handler(backend: str) -> type[ImageHandler]:
    if backend == "pillow":
        from gurobi_optimods.domino.image._pillow import PillowImageHandler

        return PillowImageHandler
    if backend == "opencv":
        from gurobi_optimods.domino.image._opencv import OpenCvImageHandler

        return OpenCvImageHandler
    if backend == "skimage":
        from gurobi_optimods.domino.image._skimage import ScikitImageHandler

        return ScikitImageHandler
    raise ValueError(
        f"'{backend}' is not a valid backend.  Please consult documentation."
    )


class TargetArray:
    """
    A class to facilitate solving the Domino Art problem by submitting a "target" array.
    Implements the Target protocol.
    """

    def __init__(self, array: npt.NDArray):
        self.array = array

    def set_params(self, **kwargs) -> None:
        warnings.warn("Ignoring parameters")

    def get_targets(self) -> npt.NDArray:
        return self.array


class TargetImage:
    """
    A class to facilitate solving the Domino Art problem by submitting a "target" array.
    Implements the Target protocol.
    """

    def __init__(self, image_handler: ImageHandler):
        self._image_handler: ImageHandler = image_handler
        self._targets: Optional[npt.NDArray] = None
        self._width: Optional[int] = None
        self._height: Optional[int] = None
        self._max_dots: Optional[int] = None

    def set_params(self, **kwargs) -> None:
        """Modifies target array according to parameters"""
        if "width" in kwargs:
            if not isinstance(kwargs["width"], int):
                raise ValueError(
                    f"Argument for 'width' parameter, {kwargs['width']}, is not an integer."
                )
            self._width = kwargs["width"]
            self._height = self._calculate_height_from_width(self._width)
        if "max_dots" in kwargs:
            max_dots = kwargs["max_dots"]
            if max_dots > 9:
                raise ValueError("The max_dots argument cannot be greater than 9.")
            self._max_dots = max_dots
        self._targets = None

    def get_targets(self) -> npt.NDArray:
        if self._targets is None:
            self._targets = self._make_targets()
        return self._targets

    def _assert_params_set(self) -> None:
        """width and max_dots must be set before we can generate targets"""
        params_missing = []
        if self._width is None:
            params_missing.append("width")
        if self._max_dots is None:
            params_missing.append("max_dots")
        if params_missing:
            raise RuntimeError(
                f"The following parameters require setting: {', '.join(params_missing)}"
            )

    def _calculate_height_from_width(self, width: int) -> int:
        orig_width, orig_height = self._image_handler.original_size
        # Need to make sure that height*width is an even number since dominoes
        #  are made up of 2 squares
        if width % 2:
            height = max(2, 2 * round(width * orig_height / orig_width / 2))
        else:
            height = max(1, round(width * orig_height / orig_width))
        return height

    def _make_targets(self) -> npt.NDArray:
        """
        Gets image as a numpy array, with correct dimensions, then scales and
        rounds values so that they range from 0, ..., max_dots
        """

        self._assert_params_set()

        assert self._width is not None
        assert self._height is not None
        # img_as_np_array = resized grayscale image (width, height)
        img_as_np_array = self._image_handler.to_numpy(self._width, self._height)

        min_, max_ = img_as_np_array.min(), img_as_np_array.max()
        if min_ == max_:
            targets = np.zeros_like(img_as_np_array)
        else:
            targets = np.round(
                self._max_dots * (max_ - img_as_np_array) / (max_ - min_)
            )

        return targets
