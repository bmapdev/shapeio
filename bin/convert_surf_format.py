#!/usr/bin/env python

"""This program converts one surface format to another. Supported formats are vtp, ucf, pial"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Ahmanson Lovelace Brain Mapping Center" \
                "University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

from shapeio import convert
import argparse


def main():
    parser = argparse.ArgumentParser(description='Convert Surface formats.')
    parser.add_argument('-surfin', dest='surfin', help='input surface [ucf,vtp]', required=True)
    parser.add_argument('-surfout', dest='surfout', help='output surface [ucf,vtp]', required=True)
    parser.add_argument('-overlay', dest='overlay', default="", help='scalar overlay [mgh]', required=False)
    parser.add_argument('-conn', dest='surfconn', default="", help='template surf with connectivity [vtp,pial]',
                        required=False)
    parser.add_argument('-applyxfm', dest='xfm', default="", help='transform file [Talairach xfm]', required=False)
    args = parser.parse_args()
    convert.surface_format(args.surfin, args.surfout, args.overlay, args.surfconn, args.xfm)
    return None

if __name__ == '__main__':
    main()

