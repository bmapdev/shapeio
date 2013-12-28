#!/usr/local/epd/bin/python

"""This program converts one surface format to another. Supported formats are vtp, ucf, pial"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Ahmanson Lovelace Brain Mapping Center" \
                "University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

import os
import sys
from shapeio import curveio
from shapeio import convert
import argparse


def main():

    parser = argparse.ArgumentParser(description='Convert Surface formats.')
    parser.add_argument('-curvein',dest='curvein',help='input curve [ucf,vtp,svg]',required=True)
    parser.add_argument('-curveout',dest='curveout',help='output curve [ucf]',required=True)
    args = parser.parse_args()
    convert.curve_format(args.curvein, args.curveout)

    return None

if __name__ == '__main__':
    main()

