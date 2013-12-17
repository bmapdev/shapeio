""" This module implements helper functions for reading and writing of surface file formats
"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi Ahmanson-Lovelace Brain Mapping Center, \
                   University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from shapeio import convert


surfin = 'shapeio/test/data/test.m'
surfin = '/Volumes/Users/sjoshi/Desktop/DE_132_02_aseg_hippo_Left_edit_smooth_1.0.nii_std_60k_jfeature.m'
surfout = 'shapeio/test/data/test.ucf'


def test_convert_surface():
    convert.surface_format(surfin, surfout)


if __name__ == "__main__":
    test_convert_surface()
