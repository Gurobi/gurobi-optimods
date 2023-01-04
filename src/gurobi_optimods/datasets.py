"""
Module for loading datasets for use in optimods examples, in the same vein
as sklearn.datasets.
"""

import pathlib

import pandas as pd

DATA_FILE_DIR = pathlib.Path(__file__).parent / "data"


class AttrDict(dict):
    """Even simpler version of sklearn's Bunch"""

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key)


def load_workforce():
    return AttrDict(
        availability=pd.read_csv(
            DATA_FILE_DIR / "availability.csv", parse_dates=["Shift"]
        ),
        pay_rates=pd.read_csv(DATA_FILE_DIR / "pay_rates.csv"),
        shift_requirements=pd.read_csv(
            DATA_FILE_DIR / "shift_requirements.csv", parse_dates=["Shift"]
        ),
    )


def load_diet():
    return AttrDict(
        categories=pd.read_csv(DATA_FILE_DIR / "diet-categories.csv"),
        foods=pd.read_csv(DATA_FILE_DIR / "diet-foods.csv"),
        nutrition_values=pd.read_csv(DATA_FILE_DIR / "diet-values.csv"),
    )
