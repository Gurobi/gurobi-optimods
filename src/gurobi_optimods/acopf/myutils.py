import sys


def breakexit(foo):
    stuff = input("("+foo+") break> ")
    if stuff == 'x' or stuff == 'q':
        sys.exit("bye")


