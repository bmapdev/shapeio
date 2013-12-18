""" This module implements file reading and writing functions for the dfs format for BrainSuite
    Also see http://brainsuite.bmap.ucla.edu for the software
"""

__author__ = "Shantanu H . Joshi"
__copyright__ = "Copyright 2013, Brandon Ayers, Ahmanson-Lovelace Brain Mapping Center, \
                 University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"


import numpy as np
import struct
import os
import sys


def readdfc(filename):

    class hdr:
        pass
    class NFV:
        pass
    Curves = []

    fid = open(filename,'rb')

    hdr.magic = np.array(struct.unpack('c'*8,fid.read(8)),dtype = 'S1')
    hdr.version = np.array(struct.unpack('c'*4,fid.read(4)),dtype = 'S1')
    hdr.hdrsize = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32') # size(int32) = 4 bytes
    hdr.dataStart = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32') # size(int32) = 4 bytes
    hdr.mdoffset = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.pdoffset  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.nContours = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')

    fid.seek(hdr.mdoffset,os.SEEK_SET)
    Mdata = struct.unpack('c'*(hdr.dataStart-hdr.mdoffset),fid.read(hdr.dataStart-hdr.mdoffset))
    hdr.xmlstr = "".join(Mdata)

    fid.seek(hdr.dataStart,os.SEEK_SET)
    for ctno in np.arange(0,hdr.nContours):
        nopts  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')

        XYZ  = np.array(struct.unpack('f'*3*nopts,fid.read(4*3*nopts)),dtype = 'float')
        XYZ = np.transpose(np.reshape(XYZ,(3,nopts),order='F'))
        Curves.append(XYZ)

    fid.close()
    if len(Curves) == 28:
        sys.stdout.write('This file ' + filename + ' is traced using a 28 curve protocol\n')
    elif len(Curves) == 26:
        sys.stdout.write('This file ' + filename + ' is traced using a 26 curve protocol\n')
    else:
        sys.stdout.write('This file ' + filename + ' is traced using an unknown protocol with ' + len(Curves) + ' number of curves. Now Exiting...\n')
        return

    return Curves, hdr