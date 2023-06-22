from typing import Tuple

import numpy.typing as npt

from gurobi_optimods.domino.image._typing import ImageHandler

try:
    import cv2
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "The opencv-python library is not installed.  Please install or change backend for image handlings."
    ) from None


class OpenCvImageHandler:
    """
    A class used to handle images with OpenCV and provide functionality to TargetImage class.
    Implements the ImageHandler protocol.
    """

    def __init__(self, img: npt.NDArray):
        """
        Parameters
        ----------
        img : npt.NDArray
            A numpy array, with 2 or 3 dimensions, as used by OpenCV for capturing
            image data

        Raises
        -------
        ValueError
            Raised when image array has neither 2 or 3 dimensions
        """
        self._img: npt.NDArray = img
        self._grayscale_img: npt.NDArray = self._make_grayscale(img)

    @staticmethod
    def _make_grayscale(img: npt.NDArray) -> npt.NDArray:
        if img.ndim == 2:
            return img
        elif img.ndim == 3:
            return cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        else:
            raise ValueError(
                f"img argument has dimension {img.ndim}.  Expected 2 or 3."
            )

    @property
    def original_size(self) -> Tuple[int, int]:
        """Returns original values of (width, height) of grayscale image"""
        shape = self._grayscale_img.shape
        return (shape[1], shape[0])

    def to_numpy(self, width: int, height: int) -> npt.NDArray:
        """Rescales original grayscale image and returns image as numpy array"""
        img = cv2.resize(self._grayscale_img, (width, height))
        return img.astype("float")

    @classmethod
    def from_path(cls, path_to_image: str) -> ImageHandler:
        """Reads image from filepath and returns instance of OpenCvImageHandler"""
        img = cv2.imread(path_to_image, cv2.IMREAD_UNCHANGED)
        return cls(img)
