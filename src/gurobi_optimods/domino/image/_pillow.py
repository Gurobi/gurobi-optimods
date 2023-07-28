from typing import Tuple

import numpy as np
import PIL

from gurobi_optimods.domino.image._typing import ImageHandler

try:
    from PIL import Image
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "The pillow library is not installed.  Please install or change backend for image handlings."
    ) from None


class PillowImageHandler:
    """
    A class used to handle images with Pillow and provide functionality to TargetImage class.
    Implements the ImageHandler protocol.
    """

    def __init__(self, img: Image):
        """
        Parameters
        ----------
        img : PIL.Image
            A color, or grayscale image
        """
        self._img: PIL.Image.Image = img
        self._grayscale_img: PIL.Image.Image = self._img.convert("L")

    @property
    def original_size(self) -> Tuple[int, int]:
        """Returns original values of (width, height) of grayscale image"""
        return self._grayscale_img.size

    # @property
    # def img(self) -> np.ndarray:
    #     return self._img

    def to_numpy(self, width, height):
        """Rescales original grayscale image and returns image as numpy array"""
        img = self._grayscale_img.resize((width, height))
        return np.asarray(img).astype("float")

    @classmethod
    def from_path(cls, path_to_image: str) -> ImageHandler:
        """Reads image from filepath and returns instance of PillowImageHandler"""
        img = Image.open(path_to_image)
        return cls(img)
