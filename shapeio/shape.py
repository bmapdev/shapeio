""" This module implements a generic shape class
"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi Ahmanson-Lovelace Brain Mapping Center, \
                   University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

import os
import surfio
import curveio
import numpy as np
import sys
import pandas


class Shape(object):

    def __init__(self):
        pass

    datatype = {'.dfs': 'surface',
                '.dfc': 'curve,',
                '.ucf': 'curve',
                '.vtp': 'surface'}

    @staticmethod
    def readfile(filename):
        filenamewoext, ext = os.path.splitext(filename)
        if ext not in Shape.datatype.keys():
            sys.stdout.write('Error: Unsupported data type. Supported data types are: ' + ', '.join(Shape.datatype.keys()))
            return

        if Shape.datatype[ext] == 'curve':
            shapobj = curveio.Curve.readfile(filename)
            return shapobj
        elif Shape.datatype[ext] == 'surface':
            shapobj = surfio.Surface.readfile(filename)
            return shapobj
        else:
            sys.stdout.write('Error: Unsupported data type. Supported data types are: ' + ', '.join(Shape.datatype.keys()))

    @staticmethod
    def writefile(filename, shapeobj):
        filenamewoext, ext = os.path.splitext(filename)
        if ext not in Shape.datatype.keys():
            sys.stdout.write('Error: Unsupported data type. Supported data types are: ' + ', '.join(Shape.datatype.keys()))
            return

        if Shape.datatype[ext] == 'curve':
            curveio.writecurve(filename, shapeobj.coords, shapeobj.attributes, shapeobj.isMultilevelUCF)
        elif Shape.datatype[ext] == 'surface':
            surfio.writesurface_new(filename, shapeobj.coords, shapeobj.faces, shapeobj.attributes, shapeobj.isMultilevelUCF)

    @staticmethod
    def read_aggregated_attributes_from_shapefilelist(shapefilelist):

        # Read first file
        s1 = Shape.readfile(shapefilelist[0])
        attributes_new = s1.attributes

        num_files = len(shapefilelist)
        attrib_size = len(attributes_new)

        attribute1_array = np.empty((num_files, attrib_size), 'float')
        attribute1_array[0, :] = attributes_new

        average_coords = s1.coords

        for i in range(1, len(shapefilelist)):
            s1 = Shape.readfile(shapefilelist[i])
            average_coords += s1.coords
            if len(s1.attributes) != attrib_size:
                sys.stdout.write("Length of attributes in Files " + shapefilelist[i] + " and " + shapefilelist[0]
                                 + " do not match. Quitting.\n")
                attribute1_array = []
                return attribute1_array
            else:
                attribute1_array[i, :] = s1.attributes

        average_coords /= len(shapefilelist)
        s1_average = s1
        s1_average.coords = average_coords

        return s1, s1_average, attribute1_array


    @staticmethod
    def read_aggregated_attributes_from_surfaces(filename):
        data_list = pandas.read_table(filename, sep='\t')
        shapefile_list = data_list['File']
        return Shape.read_aggregated_attributes_from_shapefilelist(shapefile_list)
