""" This module implements file reading and writing functions for curve formats such as ucf,
    MincObj (http://www.bic.mni.mcgill.ca)
"""
__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi Ahmanson-Lovelace Brain Mapping Center, \
                 University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

import os
import sys
import numpy as np
import dfcio
import parser
import bs4


def ReadUCFMultipleLevelsWithData(filename):
    """Read UCF file.
    Parameters
    ----------
    filename  : UCF file
    Returns
    -------
    X         : array of surfaces/curves
    atributes : float array of data values
    """
    #sys.stdout.write('Reading UCF file ' + filename+'...')
    fid = open(filename,'rt')
    dummy = ''
    # Read all the preamble till you reach contour_data
    while dummy != '<levels>':
        dummy = fid.readline().rstrip('\n')
        dummy = dummy.rstrip()
    # The next line is the number of levels
    num_levels = int(fid.readline().rstrip('\n'))

    X = []
    attributes = []
    for ctr in np.arange(0,num_levels):
        while dummy != '<point_num=>':
            dummy = fid.readline().rstrip('\n');
            dummy=dummy.rstrip()
        N = int(fid.readline().rstrip('\n'))
        dummy = fid.readline()
        Xtemp = [fid.readline().rstrip().split(' ') for i in np.arange(0,N)]

        Xtemp = np.array(Xtemp,dtype=float)

        X.append(np.array(Xtemp[:,0:3]))
        if Xtemp.shape[1] == 4:
            attributes.append(np.array(Xtemp[:,3]))

    sys.stdout.write('Done.\n')
    fid.close()
    return X,attributes


def ReadMincObjCurve(filename):

    X = []
    fid = open(filename,'rt')
    shape_type = fid.readline().split()
    if shape_type[0] == 'P':
        sys.stdout.write('This curve file ' + filename + ' contains an obj surface. Exiting without reading...\n')
        X = []
        return X
    elif shape_type[0] == 'L':
        num_pts = int(shape_type[-1])
        # Read points
        # coords = []
        coords = np.empty((3,num_pts))
        for i in range(0,num_pts):
            line_split = fid.readline().split()
            coords[:,i] = [float(line_split[0]),float(line_split[1]),float(line_split[2])]

        num_lines = int(fid.readline())
        dummy = fid.readline()
        # Read the indices of the points
        lines = [line.rstrip('\n') for line in fid.readlines()]
        lines = " ".join(lines)
        idx = [int(x) for x in lines.split()]
        idx_start = idx[0:num_lines]
        idx = idx[num_lines:]

        # Format the curves
        idx_start.insert(0,0)
        for i in range(0,num_lines):
            X.append(coords[:,range(idx_start[i],idx_start[i+1])])

        attributes = []

        return X,attributes


def readcurve(filename):

    # X = []
    # return X

    def dfc(filename):
        nCurves,hdr = dfcio.readdfc(filename)
        attributes = []
        isMultilevelUCF = True
        return nCurves,attributes,isMultilevelUCF

    def vtp(filename):
        print "to be implemented"

    def ucf(filename):
        X,attributes = ReadUCFMultipleLevelsWithData(filename)

        if len(X) > 1:
            isMultilevelUCF = True
            coords = X
            attributes = attributes
        else:
            isMultilevelUCF = False
            coords = X[0]
            if attributes:
                attributes = attributes[0]
            else:
                attributes = []


        return coords,attributes,isMultilevelUCF

    def SVG(filename):
        return ReadSVG(filename)

    def Mincobj(filename):
        return ReadMincObjCurve(filename)


    path_filename,ext = os.path.splitext(filename)
    options = {'.vtp'   : vtp,
               '.dfc'  : dfc,
               '.ucf'   : ucf,
               '.obj'   : Mincobj,
               '.svg'   : SVG,
               }
    isMultilevelUCF = False
    if ext in options:
        sys.stdout.write("Reading curve" + filename + "...")
        sys.stdout.flush()
        X,attributes,isMultilevelUCF = options[ext](filename)
        sys.stdout.write("Done.\n")
        sys.stdout.flush()
        return X,attributes,isMultilevelUCF
    else:
        sys.stdout.write("Input format " + ext + " not supported. Exiting without saving.\n")
        return None

def WriteUCF(coords,attriblabel,attributes,filename):
    T = coords.shape[0]
    attrib_flag = False
    if len(attributes) == T:
        attrib_flag = True
    else:
        sys.stdout.write("Mismatch in the length of attributes and the length of vertices of the mesh\n")
        return None

    sys.stdout.write('Writing ucf file ' + filename + '...')
    f = open(filename,'wt')
    f.write('#UCF created by surfio\n')
    f.write('<width=>\n')
    f.write('64\n')
    f.write('<height=>\n')
    f.write('146\n')
    f.write('<xrange=>\n')
    f.write("{0:f} {1:f}\n".format(0,0))
    f.write('<yrange=>\n')
    f.write("{0:f} {1:f}\n".format(0,0))
    f.write('<zrange=>\n')
    f.write("{0:f} {1:f}\n".format(0,0))
    f.write('<levels>\n');
    f.write('1\n');
    f.write('<level number=>\n');
    f.write('0.0000\n');
    f.write('<point_num=>\n');
    f.write("{0:d}\n".format(T));
    f.write('<contour_data=>\n');

    if attrib_flag:
        for i in np.arange(0,T):
            f.write("{0:f} {1:f} {2:f} {3:f}\n".format(float(coords[i,0]),float(coords[i,1]),float(coords[i,2]),float(attributes[i])))
    else:
        for i in np.arange(0,T):
            f.write("{0:f} {1:f} {2:f}\n".format(float(coords[i,0]),float(coords[i,1]),float(coords[i,2])))

    f.write('<end of level>\n')
    f.write('<end>\n');
    f.close()
    sys.stdout.write("Done.\n")

def WriteUCFMultipleLevelsWithData(filename,coords,attributes):
    # coords and attributes are now lists

    N = len(coords)
    T = len(coords[0])
    attrib_flag = False
    if len(attributes):
        L = len(attributes[0])
        if L == T:
            attrib_flag = True
        else:
            sys.stdout.write("Mismatch in the length of attributes and the length of vertices of the mesh\n")
            return None
    else:
        L = 0

    # if data is less than 3 dimensional, convert it to 3D
    #TBD

    sys.stdout.write('Writing ucf file ' + filename + '...')
    f = open(filename,'wt')
    f.write('#UCF created by surfio\n')
    f.write('<width=>\n')
    f.write('64\n')
    f.write('<height=>\n')
    f.write('146\n')
    f.write('<xrange=>\n')
    f.write("{0:f} {1:f}\n".format(0,0))
    f.write('<yrange=>\n')
    f.write("{0:f} {1:f}\n".format(0,0))
    f.write('<zrange=>\n')
    f.write("{0:f} {1:f}\n".format(0,0))
    f.write('<levels>\n');
    f.write('{0:d}\n'.format(N));

    for ii in range(0,N):
        f.write('<level number=>\n');
        f.write('{0:f}\n'.format(ii));
        f.write('<point_num=>\n');
        f.write("{0:d}\n".format(len(coords[ii])));
        f.write('<contour_data=>\n');

        if attrib_flag:
            for jj in range(0,len(coords[ii])):
                f.write("{0:f} {1:f} {2:f} {3:f}\n".format(float(coords[ii][jj][0]),float(coords[ii][jj][1]),float(coords[ii][jj][2]),float(attributes[ii][jj])))
        else:
            for jj in np.arange(0,len(coords[ii])):
                f.write("{0:f} {1:f} {2:f}\n".format(float(coords[ii][jj][0]),float(coords[ii][jj][1]),float(coords[ii][jj][2])))

        f.write('<end of level>\n')


    f.write('<end>\n');
    f.close()
    sys.stdout.write("Done.\n")
    return


def ReadSVG(filename):

    soup = bs4.BeautifulSoup(open(filename))

    # Select the maximum length curve
    all_paths = soup.findAll('path')
    maxlength = len(str(all_paths[0]))
    maxidx = 0
    for i in range(1,len(all_paths)):
        if maxlength < len(str(all_paths[i])):
            maxlength = len(str(all_paths[i]))
            maxidx = i
    maxpath = all_paths[maxidx]

    path_d = maxpath['d']

    path1 = parser.parse_path(path_d)
    X = np.empty((100,3))
    ctr = 0
    for i in np.linspace(0,1,100):
        X[ctr,0] = path1.point(i).real
        X[ctr,1] = path1.point(i).imag
        X[ctr,2] = 0.0
        ctr += 1

    attributes = []

    return X,attributes


def writecurve(filename, coords, attributes=[], isMultilevelUCF=False):

    def ucf(filename):
        if isMultilevelUCF:
            WriteUCFMultipleLevelsWithData(filename,coords,attributes)
        else:
            WriteUCF(coords,'',attributes,filename)

    path_filename,ext = os.path.splitext(filename)

    options = {'.ucf': ucf,
              }

    if ext in options:
#        sys.stdout.write("Writing surface " + filename + "...")
        options[ext](filename)
#        sys.stdout.write("Done.\n")
    else:
        sys.stdout.write("Output format " + ext + " not supported. Exiting without saving.\n")
        return None


class Curve(object):

    def __init__(self, coords=[], isMultilevelUCF=False, attributes=[]):
        self.coords = coords
        self.attributes = attributes
        self.isMultilevelUCF = isMultilevelUCF

    def read(self, filename):
        self.coords, self.attributes, self.isMultilevelUCF = readcurve(filename)

    def write(self, filename):
        writecurve(filename, self.coords, self.attributes, self.isMultilevelUCF)

    @staticmethod
    def readfile(filename):
        coords, attributes, isMultilevelUCF = readcurve(filename)
        return Curve(coords, isMultilevelUCF, attributes)
