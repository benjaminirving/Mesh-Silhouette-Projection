# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 12:36:18 2011

@author: benjamin

A mesh loads that either loads a mesh from an hdf5 file 
or loads vertices and faces
"""

import h5py
from numpy import cos, sin, zeros, sqrt
from numpy import array, dot
import numpy as np
from mesh_silhouette_cpp import *

class Mesh:
    """
    Imports a mesh and performs a number of simple operations on the mesh
    such as rotating and centering
    """
    
    def __init__(self, faces=None, vertices=None, h5file=None, 
                 h5faces=None, h5vertices=None, style="matlab"):
        
        #vertices and faces
        if (faces is not None) and (vertices is not None):
            self.vertices=np.array(vertices)
            self.faces=np.array(faces)
        
        #Loading from file (requires two database files called:
        # /vertices and /faces or the file name to be specified)
        elif h5file is not None:
            with h5py.File(h5file, 'r') as f:
                if h5vertices is None:
                    #default value to load
                    vertices=f['/vertices'].value
                else:
                    #input name
                    vertices=f['/' + h5vertices].value
                    
                if h5faces is None:
                    #default value to load
                    faces=f['/faces'].value
                else:
                    #input name
                    faces=f['/' + h5faces].value
            
            
            self.vertices=np.array(self.vertices)
            self.faces=np.array(self.faces)
            self.vertices=vertices.transpose()
            self.faces=faces.transpose()            
        else:
            print "Incorrect input!"
        
        if style is "matlab":
            #converting between matlab and python/c++
            #due to counting from 0
            self.faces=self.faces-1
            print "Changing counter from matlab to python"
        
        #Saving originals
        self.verticeso=self.vertices
        self.faceso=self.faces
                        
        #keeping track of all the rotations that are undergone (For PCA)
        #self.rot_mat_totxy = np.eye(2)
        #self.rot_mat_totxz = np.eye(2)
    
    def centre_mesh(self):
        #MESH centering
        mean1=self.vertices.sum(axis=0)/self.vertices.shape[0]
        self.vertices=self.vertices-mean1
    
    def flip_mesh(self):
        #z  = alpha + z_min
        #z' = z_max - alpha
        vz=self.vertices[:,2]
        
        vdiff = vz - vz.min()
        self.vertices[:,2] = vz.max()-vdiff
        
        if hasattr(self, 'points'):
            self.points[:,2]=vz.max() - (self.points[:,2]-vz.min())
            
        if hasattr(self, 'skel'):
            self.skel[:,2]=vz.max() - (self.skel[:,2]-vz.min())
            
    def rotate_mesh(self, angle1=1.5708, type1="xy"):
        """     Rotate mesh in x,y plane around the CM """
        # points -  vertices of mesh
        # angle1 -  angle of rotation
        #Mesh rotation
        if type1=="xy":
            ax1=(0,1)
        elif type1=="zx":
            ax1=(0,2)
        
        R=np.array([[cos(angle1), -sin(angle1)], [sin(angle1), cos(angle1)]])
        self.vertices=self.rot1(self.vertices, ax1, R)
        
        if hasattr(self, 'points'):
            self.points=self.rot1(self.points, ax1, R)       
            
        if hasattr(self, 'skel'):
            self.skel=self.rot1(self.skel, ax1, R)    
        
        
    @staticmethod
    def rot1(array1, ax1, R):
        #print array1
        array2=array1[:,ax1].transpose()
        array2=np.dot(R, array2)
        array1[:,ax1]=array2.transpose()
        return array1
        
    def normalize_mesh(self):
        frob=(self.vertices**2).sum(axis=0)
        frob=frob.sum()
        print frob
        
    def scale_mesh(self, scale_factor=1):
        self.vertices=self.vertices * scale_factor
        
    def reset_mesh(self):
        """
        Reset mesh back to the original 
        """
        self.vertices=self.verticeso
        self.faces=self.faceso

class SilhouetteMesh(Mesh):
    
    def __init__(self, *args, **kwargs):
        Mesh.__init__(self, *args, **kwargs)
       
        #silhouette subclass
        self.MshSil=MeshProject(DoubleVector(self.faces[:,0]), DoubleVector(self.faces[:,1]), 
                   DoubleVector(self.faces[:,2]), DoubleVector(self.vertices[:,0]),
                    DoubleVector(self.vertices[:,1]), DoubleVector(self.vertices[:,2]))
                    
        ue=0
                    
    def unique_edges(self):
        ue=1
        #slow bit
        print "Identification of unique edges -- slow"
        Er1=self.MshSil.unique_edges()
        print "...Done"
        
    def silhouette_vertices(self, proj_type, ViewP):
        
        try: 
        
            #updating because used to perform a rotate on the mesh
            self.MshSil.update_vertices(DoubleVector(self.vertices[:,0]), 
                                        DoubleVector(self.vertices[:,1]), 
                                        DoubleVector(self.vertices[:,2]))
                     
            #normals to the face                   
            Er2=self.MshSil.normals_face()
            
            # Finding the silhouette edges
            if proj_type == "Lodox":
                self.MshSil.silhouette_edges_lnoe(DoubleVector(ViewP))
            elif proj_type == "Normal":
                self.MshSil.silhouette_edges(DoubleVector(ViewP))
            else:
                print "Error: must choose projection type"
                
        
            #Silhouette number
            SList=self.MshSil.edge_silhouette_num()
            
            #Slist output
            F1=IntVector(SList); F2=IntVector(SList); V1=IntVector(SList); V2=IntVector(SList)
            self.MshSil.silhouette_out(F1, F2, V1, V2)
            #Faces numbper
            self.F1=np.array(F1); self.F2=np.array(F2); 
            
            #Vertices number
            self.V1=np.array(V1); self.V2=np.array(V2)
        
        except:
            print "Have unique_vertices been run?"
        
    @property
    def silvert1(self):
        self._silvert1=self.vertices[self.V1, :]
        return self._silvert1
    
    @property
    def silvert2(self):
        self._silvert2=self.vertices[self.V2, :]
        return self._silvert2
            
            
class PcaMesh(SilhouetteMesh):
    
    """
    This class allows PCA operations on a mesh.
    Derived from the basic operations on the Mesh class.
    """
    
    def input_pca(self, P, V, points=[0, 0, 0], skel=[0, 0, 0]):
        """
        P and V seem to have been flipped during import. Corrected. 
        Therefore: 
            P(m x n) m dimensions and n eigenvectors
            V(n x 1) n eigenvectors
        """

        #input P (modes x points)
        self.P=P
        self.V=V
        self.points=points
        self.points_o=self.points
#        self.points[:,[0,1]]=self.points[:,[1,0]]
        self.mnum=self.P.shape[1] # number of modes
        
        self.skel=skel
        self.skel_o=self.skel
        
    #transforming the mean mesh along one mode by amount by a
    def prin_comp_trans(self, a):
        # check existence of attributes
        # transform
        #x= x_ave + P*(mode)*amount
        
        b=zeros((self.mnum, 1))
        #print a.shape[0]
        b[0:a.shape[0]]=a
        #b=b*3*sqrt(self.V[0:25,0])
        shift=dot(self.P, b)
        
        self.vertices=self.vertices + self.three_col(shift)
        
    def scale_mesh(self, scale_factor=1):
        Mesh.scale_mesh(self, scale_factor)
        
        self.points=self.points * scale_factor
        self.skel=self.skel * scale_factor
            
    def reset_mesh(self):
        Mesh.reset_mesh(self)
        
        self.points=self.points_o
        self.skel=self.skel_o
        
        
    @staticmethod
    def three_col(x):
        #assuming single dimension
        l=x.shape[0]/3
        
        v=array([x[0:l,0], x[l:2*l,0], x[2*l:3*l,0]])
        v=v.transpose()
        
        return v
        
        #convert a single column into vector space
        
