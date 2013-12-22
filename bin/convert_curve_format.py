#!/usr/local/epd/bin/python

"""This program converts one surface format to another. Supported formats are vtp, ucf, pial"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Ahmanson Lovelace Brain Mapping Center" \
                "University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

import os
import sys
import curveio
import argparse
import numpy as np

def main():

    parser = argparse.ArgumentParser(description='Convert Surface formats.')
    parser.add_argument('-curvein',dest='curvein',help='input curve [ucf,vtp,svg]',required=True)
    parser.add_argument('-curveout',dest='curveout',help='output curve [ucf]',required=True)
    args = parser.parse_args()

    ConvertCurveFormat(args.curvein,args.curveout)

    return None

def ConvertCurveFormat(curvein,curveout):

    sin_name,sin_ext = os.path.splitext(curvein)
    sout_name,sout_ext = os.path.splitext(curveout)

    if curvein == curveout:
        sys.stdout.write("Error: Same input and output formats " + sin_ext + ". Exiting without saving anything...")
        return None

    X,attributes,isMultilevelUCF = curveio.readcurve(curvein)

    curveio.writecurve(curveout,X,attributes,isMultilevelUCF)

    return None


if __name__ == '__main__':
    main()

