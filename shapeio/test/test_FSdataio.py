""" This module implements helper functions for testing FSdataio
"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi Ahmanson-Lovelace Brain Mapping Center, \
                   University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from shapeio import FSdataio

filename = 'shapeio/test/data/bert.lh.thickness.mgh'

def test_load_mgh():
    data = FSdataio.load_mgh(filename)

if __name__ == "__main__":
    test_load_mgh()

