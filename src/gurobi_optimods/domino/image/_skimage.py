from typing import Tuple

import numpy.typing as npt

from gurobi_optimods.domino.image._typing import ImageHandler

try:
    import skimage
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "The scikit-image library is not installed.  Please install or change backend for image handlings."
    ) from None


class ScikitImageHandler:
    """
    A class used to handle images with Scikit-image and provide functionality to TargetImage class.
    Implements the ImageHandler protocol.
    """

    def __init__(self, img: npt.NDArray):
        """
        Parameters
        ----------
        img : npt.NDArray
            A numpy array, with 2 or 3 dimensions, as used by OpenCV for capturing
            image data
        """
        self._img: npt.NDArray = img
        self._grayscale_img: npt.NDArray = self._make_grayscale(img)

    @staticmethod
    def _make_grayscale(img: npt.NDArray) -> npt.NDArray:
        if img.ndim == 2:
            return img
        elif img.ndim == 3:
            return skimage.color.rgb2gray(img)
        else:
            raise ValueError(
                f"img argument has dimension {img.ndim}.  Expected 2 or 3."
            )

    @property
    def original_size(self) -> Tuple[int, int]:
        """Returns original values of (width, height) of grayscale image"""
        shape = self._grayscale_img.shape
        return (shape[1], shape[0])

    def to_numpy(self, width, height):
        """Rescales original grayscale image and returns image as numpy array"""
        img = skimage.transform.resize(
            self._grayscale_img, (height, width), anti_aliasing=True
        )
        return img.astype("float")

    @classmethod
    def from_path(cls, path_to_image: str) -> ImageHandler:
        """Reads image from filepath and returns instance of ScikitImageHandler"""
        img = skimage.io.imread(path_to_image)
        return cls(img)
