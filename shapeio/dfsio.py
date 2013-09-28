""" This module implements file reading and writing functions for the dfs format for BrainSuite
    Also see http://brainsuite.bmap.ucla.edu for the software
"""

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Brandon Ayers, Ahmanson-Lovelace Brain Mapping Center, \
                 University of California Los Angeles"
__email__ = "ayersb@ucla.edu"


import numpy as np
import struct

def readdfs(fname):
    class hdr:
        pass
    class NFV:
        pass
    fid = open(fname,'rb')
    hdr.ftype_header = np.array(struct.unpack('c'*12,fid.read(12)),dtype = 'S1')
    hdr.hdrsize = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32') # size(int32) = 4 bytes
    hdr.mdoffset = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.pdoffset  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.nTriangles  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.nVertices  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.nStrips  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.stripSize  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.normals  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.uvStart  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.vcoffset  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.labelOffset  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    hdr.vertexAttributes  = np.array(struct.unpack('i',fid.read(4)),dtype = 'int32')
    fid.seek(hdr.hdrsize) 
    NFV.faces = np.array(struct.unpack(('i'*3*hdr.nTriangles), fid.read(3*hdr.nTriangles*4 )),dtype = 'int32').reshape(hdr.nTriangles,3)
    NFV.vertices = np.array(struct.unpack('f'*3*hdr.nVertices, fid.read(3*hdr.nVertices*4)), dtype = 'float32').reshape(hdr.nVertices,3)
    if (hdr.normals>0):
        #print 'reading vertex normals.'
        fid.seek(hdr.normals)
        NFV.normals = np.array(struct.unpack('f'*3*hdr.nVertices, fid.read(3*hdr.nVertices*4)),dtype = 'float32').reshape(hdr.nVertices,3)
    if (hdr.vcoffset>0):
        #print 'reading vertex colors.'
        fid.seek(hdr.vcoffset)
        NFV.vcolor = np.array(struct.unpack('f'*3*hdr.nVertices, fid.read(3 * hdr.nVertices * 4)),dtype = 'float32').reshape(hdr.nVertices,3)
    if (hdr.uvStart>0):
        #print 'reading uv coordinates.'
        fid.seek(hdr.uvStart)
        uv = np.array(struct.unpack('f'*2*hdr.nVertices, fid.read(2*hdr.nVertices*4)),dtype = 'float32').reshape(2,hdr.nVertices)
        NFV.u = uv[0]
        NFV.v = uv[1]  
    if (hdr.labelOffset>0):
        #print 'reading vertex labels.'
        fid.seek(hdr.labelOffset)
        NFV.labels = np.array(struct.unpack('I' *hdr.nVertices, fid.read(hdr.nVertices * 2)),dtype ='uint16') #Size(unint16) = 2 bytes  
    if (hdr.vertexAttributes>0):
        #print 'reading vertex attributes.'
        fid.seek(hdr.vertexAttributes)
        NFV.attributes = np.array(struct.unpack('f' *hdr.nVertices, fid.read(hdr.nVertices * 4)), dtype = 'float32')	
    NFV.name = fname
    fid.close()
    return(NFV)
    
    
def writedfs(fname,NFV):
    ftype_header = np.array(['D', 'F', 'S', '_', 'L', 'E', ' ', 'v', '2', '.', '0','\x00']) #DFS_LEv2.0\0
    hdrsize = 184
    mdoffset = 0		# Start of metadata.
    pdoffset = 0       # Start of patient data header.
    nTriangles = len(NFV.faces.flatten()) / 3
    nVertices  = len(NFV.vertices.flatten()) / 3
    nStrips = 0
    stripSize = 0
    normals = 0
    uvoffset = 0
    vcoffset = 0
    precision = 0
    labelOffset = 0
    attributes = 0
    orientation = np.matrix(np.identity(4),dtype = 'int32')
    nextarraypos = hdrsize + 12 * (nTriangles + nVertices) # Start feilds after the header
    if (hasattr(NFV,'normals')):
        #print 'has normals'
        normals = nextarraypos
        nextarraypos = nextarraypos + nVertices * 12 #12 bytes per normal vector (3 x float32)
    if (hasattr(NFV,'vcolor')):
        #'has vcolor'
        vcoffset = nextarraypos
        nextarraypos = nextarraypos + nVertices * 12 # 12 bytes per color coordinate (3 x float32)
    if (hasattr(NFV,'u') and hasattr(NFV,'v')):
        #print 'has uv'
        uvoffset = nextarraypos
        nextarraypos = nextarraypos + nVertices *  8 # 8 bytes per uv coordinate (2 x float32)
    if (hasattr(NFV,'labels')):
        #print 'has labels'
        labelOffset = nextarraypos
        nextarraypos = nextarraypos + nVertices * 2  # 4 bytes per label (int16)
    if (hasattr(NFV,'attributes')):
        #print 'has attr'
        attributes = nextarraypos
        nextarraypos = nextarraypos + nVertices * 4  #  4 bytes per attribute (float32)
    fid = open(fname,'wb')
    fid.write(np.array(ftype_header,'S1').tostring())
    fid.write(np.array(hdrsize,'int32').tostring())
    fid.write(np.array(mdoffset,'int32').tostring())
    fid.write(np.array(pdoffset,'int32').tostring())
    fid.write(np.array(nTriangles,'int32').tostring())
    fid.write(np.array(nVertices,'int32').tostring())
    fid.write(np.array(nStrips,'int32').tostring())
    fid.write(np.array(stripSize,'int32').tostring())
    fid.write(np.array(normals,'int32').tostring())
    fid.write(np.array(uvoffset,'int32').tostring())
    fid.write(np.array(vcoffset,'int32').tostring())
    fid.write(np.array(labelOffset,'int32').tostring())
    fid.write(np.array(attributes,'int32').tostring())
    fid.write(np.zeros([1,4+15*8],'uint8').tostring()) 
    fid.write(np.array((NFV.faces),'int32').tostring()) 
    fid.write(np.array(NFV.vertices,'float32').tostring()) 
    if (normals>0):
        #print 'writing normals' 
        fid.write(np.array(NFV.normals,'float32').tostring())
    if vcoffset> 0:
        #print 'writing color'
        fid.write(np.array(NFV.vColor,'float32').tostring())
    if uvoffset > 0:
        #print 'writing uv'
        fid.write(np.array([NFV.u.flatten(), NFV.v.flatten()],'float32').tostring())
    if (labelOffset > 0):
        #print 'writing labels' 
        fid.write(np.array(NFV.labels,'int16').tostring())
    if (attributes > 0):
        #print 'writing attributes' 
        fid.write(np.array(NFV.attributes,'float32').tostring())
    fid.close()
