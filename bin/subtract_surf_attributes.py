#!/usr/bin/env python

"""This program performs subtracts the attributes of surface 2 from surface 1. Supported formats are vtp, ucf, pial"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Ahmanson Lovelace Brain Mapping Center" \
                "University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from shapeio.shape import Shape
from sys import stdout
import argparse


def main():
    parser = argparse.ArgumentParser(description='Subtract surface attributes (surface 1 - surface 2). '
                                                 'Output surface attributes will contain subtracted attributes.'
                                                 'Output surface vertices will be the same as the first surface.')
    parser.add_argument('-s1', dest='s1file', help='input surface 1 [ucf,vtp]', required=True)
    parser.add_argument('-s2', dest='s2file', help='input surface 2 [ucf,vtp]', required=True)
    parser.add_argument('-o', dest='soutfile', help='output surface [ucf,vtp]', required=True)
    args = parser.parse_args()
    subtract_surf_attributes(args.s1file, args.s2file, args.soutfile)
    return None


def subtract_surf_attributes(s1file, s2file, soutfile):
    s1 = Shape.readfile(s1file)
    s2 = Shape.readfile(s2file)

    if s1.attributes.any() and s2.attributes.any():
        if len(s1.attributes) != len(s2.attributes):
            stdout.write('Error: Attribute lengths do not match. Quitting without subtracting.\n')
            return
        else:
            s1.attributes = s1.attributes - s2.attributes
            Shape.writefile(soutfile, s1)
    else:
        stdout.write('Error: Surfaces do not contain attributes. Quitting without subtracting.\n')
        return

if __name__ == '__main__':
    main()
