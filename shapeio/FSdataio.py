""" This module implements readers and writers of FS data formats
    Based on Freesurfer (c) MGH file formats
"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi Ahmanson-Lovelace Brain Mapping Center, \
                   University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

import numpy as np
from struct import unpack


def load_mgh(filename):

    fobj = open(filename, "rb")
    ver = unpack('>i', fobj.read(4))[0]
    dim1 = unpack('>i', fobj.read(4))[0]
    dim2 = unpack('>i', fobj.read(4))[0]
    dim3 = unpack('>i', fobj.read(4))[0]
    frames = unpack('>i', fobj.read(4))[0]
    datatype = unpack('>i', fobj.read(4))[0]

    type_info = {0: (1, 'c'),
                 1: (4, 'i'),
                 3: (4, 'f'),
                 4: (2, 'h'),
                 }
    bytes_per_voxel, fmt_string = type_info[datatype]
    fobj.seek(284)  # Skip the header and jump past the unused bytes. See load_mgh supplied by Freesurfer (c) MGH
    nvoxels = dim1 * dim2 * dim3 * frames
    data = np.array(unpack('>' + fmt_string*nvoxels, fobj.read(nvoxels*bytes_per_voxel)), dtype='float32')
    return data
