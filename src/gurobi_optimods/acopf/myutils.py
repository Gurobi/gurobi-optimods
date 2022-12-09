import sys


def breakexit(foo):
    """Will be removed in the final version"""
    stuff = input("("+foo+") break> ")
    if stuff == 'x' or stuff == 'q':
        sys.exit("bye")


