""" This module implements file reading and writing functions for the dfs format for BrainSuite
    Also see http://brainsuite.bmap.ucla.edu for the software
"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Ahmanson-Lovelace Brain Mapping Center, \
                 University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from shapeio.surfio import Surface

dfsfile = 'data/left_test.dfs'


def test_read():
    s1 = Surface.readfile(dfsfile)


def test_write():
    pass


if __name__ == "__main__":
    test_read()
    test_write()


