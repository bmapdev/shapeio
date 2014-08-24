""" This module implements file reading and writing functions for vtk and vtp surface formats (www.vtk.org)
"""

__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi Ahmanson-Lovelace Brain Mapping Center, \
                 University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"

import sys
import numpy as np
import vtk
from vtk.util import numpy_support


def ReadVTK_Polydata(vtkfile):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(vtkfile)
    reader.Update()

    mesh = reader.GetOutput()
    coords = numpy_support.vtk_to_numpy(mesh.GetPoints().GetData())
    pointdata = mesh.GetPointData()

    attributes = []
    # TODO Perhps return the array by name in the future
    if pointdata.GetNumberOfArrays() >= 1:  # Attributes present
        attributes = numpy_support.vtk_to_numpy(pointdata.GetArray(0))

    temp_faces = numpy_support.vtk_to_numpy(mesh.GetPolys().GetData())

    num_faces = mesh.GetPolys().GetNumberOfCells()
    strides = np.arange(0, len(temp_faces), 4)
    faces = np.delete(temp_faces, strides)
    faces = np.reshape(faces, (num_faces, 3))

    return coords, faces, attributes


def ReadVTK_XML_Polydata(vtkfile):

    sys.stdout.write("Reading VTK XML Polydata (Assuming triangular mesh)...")
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtkfile)
    reader.Update()
    mesh = reader.GetOutput()
    coords = numpy_support.vtk_to_numpy(mesh.GetPoints().GetData())
    pointdata = mesh.GetPointData()

    attributes = []
    # TODO Perhps return the array by name in the future
    if pointdata.GetNumberOfArrays() >= 1:  # Attributes present
        if pointdata.HasArray('VoxelData'):  # If VoxelData present use it or else pick the first array
            attributes = numpy_support.vtk_to_numpy(pointdata.GetArray('VoxelData'))
        else:
            attributes = numpy_support.vtk_to_numpy(pointdata.GetArray(0))

    temp_faces = numpy_support.vtk_to_numpy(mesh.GetPolys().GetData())

    num_faces = mesh.GetPolys().GetNumberOfCells()
    strides = np.arange(0, len(temp_faces), 4)
    faces = np.delete(temp_faces, strides)
    faces = np.reshape(faces, (num_faces, 3))

    return coords, faces, attributes


def write_vtk_xml_polydata_curve(filename, coords, attributes = []):
    if coords.shape[0] < coords.shape[1]:
        coords = coords.T

    coords = coords.copy()  # To ensure C_CONTIGUOUS is True

    polydata = vtk.vtkPolyData()
    points = vtk.vtkPoints()

    points.SetData(numpy_support.numpy_to_vtk(coords))
    polydata.SetPoints(points)

    lines = vtk.vtkCellArray()
    lines.InsertNextCell(coords.shape[0])
    for i in np.arange(coords.shape[0]):
        lines.InsertCellPoint(i)

    polydata.SetLines(lines)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInput(polydata)
    writer.SetFileName(filename)
    writer.SetDataModeToAscii()
    writer.Write()


def write_vtk_xml_polydata_curve_set(filename, coords_set, attributes = []):
    levels = len(coords_set)


    vtkmultipieceobj = vtk.vtkMultiPieceDataSet()
    vtkmultipieceobj.SetNumberOfPieces(levels)

    for ii in np.arange(levels):

        coords = coords_set[ii]

        if coords.shape[0] < coords.shape[1]:
            coords = coords.T

        coords = coords.copy()  # To ensure C_CONTIGUOUS is True

        polydata = vtk.vtkPolyData()
        points = vtk.vtkPoints()

        points.SetData(numpy_support.numpy_to_vtk(coords))
        polydata.SetPoints(points)

        lines = vtk.vtkCellArray()
        lines.InsertNextCell(coords.shape[0])
        for i in np.arange(coords.shape[0]):
            lines.InsertCellPoint(i)

        polydata.SetLines(lines)
        vtkmultipieceobj.SetPiece(ii, polydata)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInput(vtkmultipieceobj)
    writer.SetFileName(filename)
    writer.SetDataModeToAscii()
    writer.Write()


def write_multilevel_polyline_to_vtp(filename, coords_set, attributes = []):
    levels = len(coords_set)
    fid = open(filename, mode='wt')
    fid.write('<?xml version="1.0"?>\n')
    fid.write('<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n')
    fid.write('<PolyData>\n')


    T = coords_set[0].shape[0]
    idx = np.arange(T)
    g = [''.join(str(i)) for i in idx]
    idx_str = ' '.join(g)

    for ii in np.arange(levels):
        coords = coords_set[ii]
        fid.write('<Piece NumberOfPoints="{0}" NumberOfVerts="0" NumberOfLines="1" NumberOfStrips="0" NumberOfPolys="0">\n'.format(T))
        fid.write('<Points>\n')
        fid.write('<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        for i in coords:
            fid.write("{0:f} {1:f} {2:f}\n".format(float(i[0]), float(i[1]), float(i[2])))

        fid.write('</DataArray>\n')
        fid.write('</Points>\n')
        fid.write('<Lines>\n')
        fid.write('<DataArray type="Int32" Name="connectivity" format="ascii">\n')

        fid.write(idx_str + '\n')
        fid.write('\n</DataArray>\n')
        fid.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
        fid.write(str(T) + '\n')
        fid.write('</DataArray>\n')
        fid.write('</Lines>\n')
        fid.write('</Piece>\n')

    fid.write('</PolyData>\n')
    fid.write('</VTKFile>\n')
    fid.close()


def WriteVTK_XML_Polydata_old(vtkfile,coords,faces,attriblabel,attrib,ascii_flag=True):

    sys.stdout.write('Writing vtp file ' + vtkfile + '...')
    T = coords.shape[0]
    F = faces.shape[0]

    if len(attrib) > 0 and len(attrib) != coords.shape[0]:
        print "Attribute length " + str(len(attrib)) + " not the same as the length of coordinate vertices " + str(coords.shape[0]) +  " . Not saving file"
        return None

    fid = open(vtkfile,mode='wt')
    fid.write('<?xml version="1.0"?>\n')
    fid.write('<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n');
    fid.write('<PolyData>\n');
    fid.write('<Piece NumberOfPoints="{0}" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="{1}">\n'.format(T,F))
    fid.write('<Points>\n');
    fid.write('<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');

    for i in coords:
        fid.write("{0:f} {1:f} {2:f}\n".format(float(i[0]),float(i[1]),float(i[2])))

    fid.write('</DataArray>\n');
    fid.write('</Points>\n');

    if len(attrib) > 0:
        if attriblabel == "":
            attriblabel = "attrib-label"

        fid.write('<PointData Scalars="{0}">\n'.format(attriblabel))
        fid.write('<DataArray type="Float32" Name="{0}" format="ascii">\n'.format(attriblabel))
        for i in attrib:
            fid.write("{0:f} ".format(float(i)))

        fid.write('\n</DataArray>\n')
        fid.write('</PointData>\n')

    fid.write('<Polys>\n');
    fid.write('<DataArray type="Int32" Name="connectivity" format="ascii">\n');

    if faces.min() == 1:
        faces -= 1

    triangles = np.ravel(faces,'C').astype('int')

    for i in triangles:
        fid.write("{0:d} ".format(i))


    fid.write('\n</DataArray>\n');
    fid.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')

    idx_array = np.arange(3,3*faces.shape[0]+3,3)

    for i in idx_array:
        fid.write("{0:d} ".format(i))

    fid.write('\n</DataArray>\n')
    fid.write('</Polys>\n')

    fid.write('</Piece>\n')
    fid.write('</PolyData>\n')
    fid.write('</VTKFile>\n')

    fid.close()
    sys.stdout.write("Done.\n")


    return None
