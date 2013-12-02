""" This module implements helper functions for reading and writing of surface file formats
"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi Ahmanson-Lovelace Brain Mapping Center, \
                   University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"


import os
import sys
import numpy as np
import vtkio
import dfsio
import re
import StringIO
import pandas as pd

# Optionally import the parse module
import imp
try:
    import parse
except ImportError:
    parse_module_exists = False



FORMATS_TYPES = {
    'd': int,
    'f': float,
    's': str,
}


FORMATS_REGEXES = {
    'd': re.compile(r'(?:\s|\b)*([+-]?\d+)(?:\s|\b)*'),
    'f': re.compile(r'(?:\s|\b)*([+-]?\d+\.?\d*)(?:[eE][+-]?\d+)?(?:\s|\b)*'),
    's': re.compile(r'\b(\w+)\b'),
}


FORMAT_FIELD_REGEX = re.compile(r'%(s|d|f)')


def scan_input(format_string, stream, max_size=float('+inf'), chunk_size=1024):
    """Scan an input stream and retrieve formatted input."""

    chunk = ''
    format_fields = format_string.split()[::-1]
    while format_fields:
        fields = FORMAT_FIELD_REGEX.findall(format_fields.pop())
        if not chunk:
            chunk = _get_chunk(stream, chunk_size)

        for field in fields:
            field_regex = FORMATS_REGEXES[field]
            match = field_regex.search(chunk)
            length_before = len(chunk)
            while match is None or match.end() >= len(chunk):
                chunk += _get_chunk(stream, chunk_size)
                if not chunk or length_before == len(chunk):
                    if match is None:
                        raise ValueError('Missing fields.')
                    break
            text = match.group(1)
            yield FORMATS_TYPES[field](text)
            chunk = chunk[match.end():]



def _get_chunk(stream, chunk_size):
    try:
        return stream.read(chunk_size)
    except EOFError:
        return ''

def ReadPial(filename):
    fid = open(filename,'rb')
    b1, b2, b3 = np.fromfile(fid, dtype=np.uint8, count=3)

    magic_code = (b1 << 16) + (b2 << 8) + b3
    if magic_code == 16777214:  # Currently only read triangular meshes
        dummy = fid.readline()
        dummy = fid.readline()
        vnum = np.fromfile(fid,dtype=">i4",count=1)[0]
        fnum = np.fromfile(fid,dtype=">i4",count=1)[0]
        coords = np.fromfile(fid,dtype=">f4",count=vnum*3).reshape(vnum,3)
        faces = np.fromfile(fid,dtype=">i4",count=fnum*3).reshape(fnum,3)
    else:
        sys.stdout.write('Invalid file type. Currently only Freesurfer triangular meshes are supported.')
    return coords, faces



def WriteUCF(coords,attriblabel,attributes,filename):
    T = coords.shape[0]
    attrib_flag = False
    if len(attributes):
        L = attributes.shape[0]
        if L == T:
            attrib_flag = True
        else:
            sys.stdout.write("Mismatch in the length of attributes and the length of vertices of the mesh\n")
            return None
    else:
        L = 0

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


def ReadMincObj(filename):
    fid = open(filename,'rt')


    shape_type = fid.readline().split()

    if shape_type[0] == 'P':
        num_pts = int(shape_type[-1])

        # Read vertices
        coords = []
        for i in range(0,num_pts):
            line_split = fid.readline().split()
            Xtemp = [float(line_split[0]),float(line_split[1]),float(line_split[2])]
            coords.append(Xtemp)

        # Read a blank line
        fid.readline()

        # Read normals
        normals = []
        for i in range(0,num_pts):
            line_split = fid.readline().split()
            Xtemp = [float(line_split[0]),float(line_split[1]),float(line_split[2])]
            normals.append(Xtemp)

        # Read a blank line
        fid.readline()
        # Read number of triangles
        line = fid.readline()
        num_triangles = int(line.split()[0])
        line = fid.readline()
        color_type = int(line.split()[0])

        if color_type == 2:
        # Read colors
            pt_colors = []
            for i in range(0,num_pts):
                line_split = fid.readline().split()
                Xtemp = [float(line_split[0]),float(line_split[1]),float(line_split[2]),float(line_split[3])]
                pt_colors.append(Xtemp)
        elif color_type == 0:
            line = fid.readline()
            line_split = line.split()
            pt_colors = [float(line_split[0]),float(line_split[1]),float(line_split[2]),float(line_split[3])]

        # Read a blank line
        dummy = fid.readline()

        # Till the next blank line, read all lines
        dummy = fid.readline()
        while dummy != '':
            dummy = fid.readline().rstrip('\n')

        faces = []
        for i in range(0,num_triangles):
            line_split = fid.readline().split()
            Xtemp = [int(line_split[0]),int(line_split[1]),int(line_split[2])]
            faces.append(Xtemp)

        fid.close()
        isMultilevelUCF = False
        coords = np.array(coords)
        attributes = pt_colors
        return coords, faces, attributes, isMultilevelUCF


    elif shape_type[0] == 'L': # This means the obj file contains lines
        sys.stdout.write('This surface file ' + filename + ' contains obj curves. Exiting without reading...\n')
        coords =  faces = attributes = isMultilevelUCF = False
        return coords, faces, attributes, isMultilevelUCF

def ReadCCBBM_sphere(filename):
    fid = open(filename,'rt')
    lines = fid.readlines()
    coords = []

    # Skip all lines until vertex
    ctr = 0
    while lines[ctr][0] == '#':
        ctr += 1

    # First read the vertices
    line = lines[ctr]
    line_split = line.split()
    ctr += 1
    # ctr = 1;
    while line_split[0] == 'Vertex':
        s1 = StringIO.StringIO(lines[ctr])
        # data = surfio.scan_input('Vertex %d %f %f %f {min=(%f) coord=(%f %f %f) max=(%f) cn=(%f) gauss=(%f) mean=(%f) sphere=(%f %f %f) mean=(%f) si=(%f)}',s1)
        data = scan_input('Vertex %d %f %f %f {coord=(%f %f %f) min=(%f) max=(%f) cn=(%f) normal=(%f %f %f) gauss=(%f) mean=(%f) sphere=(%f %f %f) dist=(%f) si=(%f)}',s1);
        for d in data:
            sys.stdout.write('{0} '.format(repr(d)))
            # print repr(d)

        sys.stdout.write('\n')
        ctr = ctr + 1


        # Xtemp = [float(line_split[2]),float(line_split[3]),float(line_split[4])]
        # coords.append(Xtemp)
        # line = lines[ctr]
        # line_split = line.split()
        # ctr += 1

    coords = np.array(coords)


def readccbbm(filename):
    fid = open(filename, 'rt')
    lines = fid.readlines()
    coords = []

    # Skip all lines until vertex
    ctr = 0
    while lines[ctr][0:6] != 'Vertex':
        ctr += 1

    # First read the vertices
    line = lines[ctr]
    line_split = line.split()
    ctr += 1
    # ctr = 1;
    attributes = []
    while line_split[0] == 'Vertex':
        vtx1 = float(line_split[2])
        vtx2 = float(line_split[3])
        vtx3 = float(line_split[4])
        radial_dist = float(re.findall('\d*\.?\d+', line_split[5])[0])
        coords.append([vtx1, vtx2, vtx3])
        attributes.append(radial_dist)
        line = lines[ctr]
        line_split = line.split()
        ctr += 1

    coords = np.array(coords)
    # The rest of the lines are faces
    faces = []
    ctr -= 1
    for ii in range(ctr, len(lines)):
        line_split = lines[ii].split()
        faces.append([int(line_split[2]), int(line_split[3]), int(line_split[4])])

    faces = np.array(faces)
    if faces.min() == 1:
        faces -= 1

    isMultilevelUCF = False
    return coords, faces, attributes, isMultilevelUCF



def readccbbm_using_parse(filename):

    fid = open(filename, 'rt')
    lines = fid.readlines()
    regex=r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'

    coords = []

    # Skip all lines until vertex
    ctr = 0
    parser1 = parse.compile("Vertex {:d} {:g} {:g} {:g} {Jfeature=({:g} {:g} {:g} {:g} {:g} {:g} {:g})}")
    parser2 = parse.compile("Vertex {:d} {:g} {:g} {:g}")

    while lines[ctr][0:6] != 'Vertex':
        ctr += 1

    # First read the vertices
    line = lines[ctr]
    line_split = line.split()
    ctr += 1
    # ctr = 1;
    attributes = []
    while line_split[0] == 'Vertex':
        result = parser1.parse(line)
        if result is not None:
            (idx, vtx1, vtx2, vtx3, radial_dist, mTBM1, mTBM2, mTBM3, detJacobian, eig1Jacobian, eig2Jacobian) = result.fixed
            attributes.append(radial_dist)
        else:
            result = parser2.parse(line)
            if result is not None:
                (idx, vtx1, vtx2, vtx3) = result.fixed
            else:
                sys.stdout.write('Cannot parse the line' + line)

        coords.append([vtx1, vtx2, vtx3])
        line = lines[ctr]
        line_split = line.split()
        ctr += 1

    coords = np.array(coords)
    # The rest of the lines are faces
    faces = []
    ctr -= 1
    for ii in range(ctr, len(lines)):
        line = lines[ii]
        result = parse.search("Face {:d} {:d} {:d} {:d}", line)
        (idx, face1, face2, face3) = result.fixed
        faces.append([face1, face2, face3])

    faces = np.array(faces)
    if faces.min() == 1:
        faces -= 1

    isMultilevelUCF = False
    return coords, faces, attributes, isMultilevelUCF

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


def Read(filename):

    coords,faces,attributes = readsurface(filename)
    return Surface(coords,faces,attributes)

def readsurface_new(filename):

    def Mincobj(filename):
        return ReadMincObj(filename)

    def dfs(filename):
        NFV = dfsio.readdfs(filename)
        coords = NFV.vertices
        faces = NFV.faces
        attributes = []
        # attributes = NFV.attributes
        isMultilevelUCF = False
        return coords,faces,attributes,isMultilevelUCF

    def vtp(filename):
        coords,faces,attributes = vtkio.ReadVTK_XML_Polydata(filename)
        isMultilevelUCF = False
        return coords,faces,attributes,isMultilevelUCF

    def ccbbm(filename):
        return readccbbm(filename)

    def ucf(filename):
        X,attributes = ReadUCFMultipleLevelsWithData(filename)

        if len(X) > 1:
            isMultilevelUCF = True
            coords = X
            faces = []
            attributes = attributes
        else:
            isMultilevelUCF = False
            coords = X[0]
            faces = []
            attributes = attributes[0]

        return coords,faces,attributes,isMultilevelUCF

    def vtk(filename):
        coords,faces,attributes = vtkio.ReadVTK_Polydata(filename)
        isMultilevelUCF = False
        return coords,faces,attributes,isMultilevelUCF

    def pial(filename):
        # coords,faces = FSio.read_geometry(filename)
        coords,faces = ReadPial(filename)
        attributes = []
        isMultilevelUCF = False
        return coords,faces,attributes,isMultilevelUCF

    def FSsphere(filename):
        # coords,faces = FSio.read_geometry(filename)
        coords,faces = ReadPial(filename)
        attributes = []
        isMultilevelUCF = False
        return coords,faces,attributes,isMultilevelUCF

    path_filename,ext = os.path.splitext(filename)
    options = {'.vtp'   : vtp,
               '.ucf'   : ucf,
               '.dfs'   : dfs,
               '.vtk'   : vtk,
               '.pial'  : pial,
               '.reg'   : FSsphere,
               '.sphere': pial,
               '.m'     : ccbbm,
               '.obj'   : Mincobj,
               }

    isMultilevelUCF = False
    if ext in options:
        sys.stdout.write("Reading surface " + filename + "...")
        sys.stdout.flush()
        coords,faces,attributes,isMultilevelUCF = options[ext](filename)
        sys.stdout.write("Done.\n")
        sys.stdout.flush()
        return coords,faces,attributes,isMultilevelUCF
    else:
        sys.stdout.write("Input format " + ext + " not supported. Exiting without saving.\n")
        return None


def readsurface(filename):

    def vtp(filename):
        return vtkio.ReadVTK_XML_Polydata(filename)

    def ucf(filename):
        X,attributes = ReadUCFMultipleLevelsWithData(filename)

        if len(X) > 1:
            isMultilevelUCF = True
            coords = X
            faces = []
            attributes = attributes
        else:
            coords = X[0]
            faces = []
            attributes = attributes[0]

        return coords,faces,attributes

    def vtk(filename):
        return None

    def pial(filename):
        # coords,faces = FSio.read_geometry(filename)
        coords,faces = ReadPial(filename)
        attributes = []
        return coords,faces,attributes

    def FSsphere(filename):
        # coords,faces = FSio.read_geometry(filename)
        coords,faces = ReadPial(filename)
        attributes = []
        return coords,faces,attributes

    path_filename,ext = os.path.splitext(filename)
    options = {'.vtp'   : vtp,
               '.ucf'   : ucf,
               '.vtk'   : vtk,
               '.pial'  : pial,
               '.reg'   : FSsphere,
               '.sphere': pial,
              }

    if ext in options:
        sys.stdout.write("Reading surface " + filename + "...")
        sys.stdout.flush()
        coords,faces,attributes = options[ext](filename)
        sys.stdout.write("Done.\n")
        sys.stdout.flush()
        return coords,faces,attributes
    else:
        sys.stdout.write("Input format " + ext + " not supported. Exiting without saving.\n")
        return None

def Write(filename,surf):
    writesurface(filename,surf.coords,surf.faces,surf.attributes)



def writesurface(filename,coords,faces,attributes):

    def vtp(filename):
        vtkio.WriteVTK_XML_Polydata_old(filename,coords,faces,'',attributes)

    def ucf(filename):
        WriteUCF(coords,'',attributes,filename)

    def vtk(filename):
        sys.stdout.write("Not implemented.")
        return None

    def pial(filename):
        # FSio.write_geometry(filename,coords,faces)
        # sys.stdout.write("Not implemented.")
        return None

    path_filename,ext = os.path.splitext(filename)

    options = {'.vtp'  : vtp,
               '.ucf'  : ucf,
               '.vtk'  : vtk,
               '.pial' : pial,
               }
    if ext in options:
    #        sys.stdout.write("Writing surface " + filename + "...")
        options[ext](filename)
    #        sys.stdout.write("Done.\n")
    else:
        sys.stdout.write("Output format " + ext + " not supported. Exiting without saving.\n")
        return None

def writesurface_new(filename,coords,faces,attributes=[],isMultilevelUCF=False):

    def vtp(filename):
        vtkio.WriteVTK_XML_Polydata_old(filename,coords,faces,'',attributes)

    def dfs(filename):
        s1 = Surface()
        s1.vertices = coords
        s1.faces = faces
        # NFV.attributes = attributes
        dfsio.writedfs(filename,s1)


    def ucf(filename):
        if isMultilevelUCF:
            WriteUCFMultipleLevelsWithData(filename,coords,attributes)
        else:
            WriteUCF(coords,'',attributes,filename)

    def vtk(filename):
        sys.stdout.write("Not implemented.")
        return None

    def pial(filename):
        # FSio.write_geometry(filename,coords,faces)
        # sys.stdout.write("Not implemented.")
        return None

    path_filename,ext = os.path.splitext(filename)

    options = {'.vtp'  : vtp,
               '.ucf'  : ucf,
               '.vtk'  : vtk,
               '.pial' : pial,
    }
    if ext in options:
#        sys.stdout.write("Writing surface " + filename + "...")
        options[ext](filename)
#        sys.stdout.write("Done.\n")
    else:
        sys.stdout.write("Output format " + ext + " not supported. Exiting without saving.\n")
        return None

def convertUCFSurfacetoVTP(ucf_surf,surf_with_conn,vtp_surf):
    X,attributes = ReadUCFMultipleLevelsWithData(ucf_surf)
    coords = np.array(X[0])

    # conn_coords,conn_faces = FSio.read_geometry(surf_with_conn)
    conn_coords,conn_faces = ReadPial(surf_with_conn)
    vtkio.WriteVTK_XML_Polydata_old(vtp_surf,X[0],conn_faces,'',attributes[0])


def aggregate_ucf_attributes(attributes):

    N = len(attributes)
    T = len(attributes[0])

    attribues_aggregate = []
    for i in range(0,N):
        attribues_aggregate = np.append(attribues_aggregate,attributes[i])

    return attribues_aggregate


def read_aggregated_attributes_from_surfaces(filename):
    data_list = pd.read_table(filename, sep='\t')

    # Read first file
    s1 = Surface()
    s1.read(data_list['File'][0])
    if s1.ismultilevelUCF:
        attributes_new = aggregate_ucf_attributes(s1.attributes)
    else:
        attributes_new = s1.attributes

    num_files = len(data_list['File'])
    attrib_size = len(attributes_new)

    attribute1_array = np.empty((num_files, attrib_size), 'float')
    attribute1_array[0, :] = attributes_new

    file_list = data_list['File']
    average_coords = s1.coords

    if s1.ismultilevelUCF:
        for i in range(1, len(file_list)):
            s1.read(data_list['File'][i])
            attributes_new = aggregate_ucf_attributes(s1.attributes)
            if len(attributes_new) != attrib_size:
                sys.stdout.write("Length of attributes in Files " + data_list['File'][i] + " and " + data_list['File'][0]
                                 + " do not match. Quitting.\n")
                attribute1_array = []
                return attribute1_array
            else:
                attribute1_array[i, :] = attributes_new
    else:
        for i in range(1, len(file_list)):
            s1.read(data_list['File'][i])
            average_coords += s1.coords
            if len(s1.attributes) != attrib_size:
                sys.stdout.write("Length of attributes in Files " + data_list['File'][i] + " and " + data_list['File'][0]
                                 + " do not match. Quitting.\n")
                attribute1_array = []
                return attribute1_array
            else:
                attribute1_array[i, :] = s1.attributes

    average_coords /= len(file_list)
    s1_average = s1
    s1_average.coords = average_coords

    return s1, s1_average, attribute1_array


def write_pvalues_to_surface(pvalue_array, s1, filename):

    if s1.ismultilevelUCF is True:
        pvalue_array = np.reshape(pvalue_array, (len(s1.attributes), len(s1.attributes[0])), 'C')
        pvalues = []
        for ii in range(0, len(s1.attributes)):
            pvalues.append(np.array(pvalue_array[ii, :]))
    else:
        pvalues = pvalue_array

    s1.attributes = pvalues
    s1.write(filename)


class Surface(object):

    def __init__(self,coords=[],faces=[],ismultilevelUCF=False,attributes=[],attriblabel=[]):
        self.coords = coords
        self.faces = faces
        self.attributes = attributes
        self.attriblabel = attriblabel
        self.ismultilevelUCF = ismultilevelUCF

    def read(self,filename):
        self.coords,self.faces,self.attributes,self.ismultilevelUCF = readsurface_new(filename)

    def write(self,filename):
        writesurface_new(filename,self.coords,self.faces,self.attributes,self.ismultilevelUCF)

    @staticmethod
    def readfile(filename):
        coords,faces,attributes,ismultilevelUCF = readsurface_new(filename)
        return Surface(coords,faces,ismultilevelUCF,attributes,"")

    def computeBoundingSphere(self):

        minX,minY,minZ = np.min(self.coords[:,0]),np.min(self.coords[:,1]),np.min(self.coords[:,2])
        maxX,maxY,maxZ = np.max(self.coords[:,0]),np.max(self.coords[:,1]),np.max(self.coords[:,2])
        xc,yc,zc = np.mean(self.coords[:,0]),np.mean(self.coords[:,1]),np.mean(self.coords[:,2])
        radius = np.sqrt( (minX - maxX)**2 + (minY - maxY)**2 + (minZ - maxZ)**2 )*0.5

        return xc,yc,zc,radius

