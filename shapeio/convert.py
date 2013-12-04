""" This module implements helper functions for reading and writing of surface file formats
"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi Ahmanson-Lovelace Brain Mapping Center, \
                   University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

import os
import sys
from shapeio import surfio
from shapeio import curveio
import numpy as np


def surface_format(surfin, surfout, surf_conn="", xfm=""):
    sin_name, sin_ext = os.path.splitext(surfin)

    if surfin == surfout:
        sys.stdout.write("Error: Same input and output formats " + sin_ext + ". Exiting without saving anything...")
        return None

    in_faces = []
    if sin_ext == ".ucf":
        if surf_conn == "":
            sys.stdout.write('Attempting to convert ucf (point set) to vtp (triangular mesh). '
                             'A template connectivity surface is required. Exiting without saving...')
            return None
        else:
            conn_coords, conn_faces, conn_attributes, isMultilevelUCF = surfio.readsurface_new(surf_conn)
            in_faces = conn_faces

    in_coords, temp_faces, in_attributes, isMultilevelUCF = surfio.readsurface_new(surfin)
    if not len(in_faces):
        in_faces = temp_faces

    if xfm != "":
        sys.stdout.write("Applying transform..." + xfm + "\n")
        mtx = np.loadtxt(xfm)
        print mtx
        in_coords = np.dot(np.c_[in_coords, np.ones(len(in_coords))], mtx.T)[:, :3]
        sys.stdout.write("Done.\n")

    surfio.writesurface_new(surfout, in_coords, in_faces, in_attributes)

    return None


def curve_format(curvein, curveout):

    sin_name, sin_ext = os.path.splitext(curvein)

    if curvein == curveout:
        sys.stdout.write("Error: Same input and output formats " + sin_ext + ". Exiting without saving anything...")
        return None

    X, attributes, isMultilevelUCF = curveio.readcurve(curvein)

    curveio.writecurve(curveout, X, attributes, isMultilevelUCF)

    return None

