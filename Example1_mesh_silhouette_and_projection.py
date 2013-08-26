"""
run_pcaproj
Main file for running the projection algorithm
1) import statistical model
2) project statistical model to 2D
3) align statistical model with 2D airway
4) optimise the fit of the statistical components

Benjamin Irving
2013 / 08 / 25

"""
print __doc__

#Python modules
import numpy as np
import time, copy, os, pickle
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.ndimage as nd

#My modules
import Optimisation_Fitting as ofit
import point_man_class as pmc
import mesh3D_mod.project_module as pm
import mesh3D_mod.mesh_class as mc
from mesh3D_mod import meshplot_module
from importdata.matinteract_module import H5fileData
import image2D.imagefilt2D_module as ft
import mayavi.mlab as mlab
        
    
# Viewing direction
#Lodox source to detector distance 1299mm
# Projecting onto a plane
ViewP=[1000.0, 0.0, 0.0]
p0=[-100.0, 0.0, 0.0]
proj_type="Lodox"

# reading in mesh and statistical model
math5=H5fileData("data_model/meshdata.h5")
scale1=1.0/math5.scale1.mean()

#Full example mesh
MHexam=mc.SilhouetteMesh(faces=math5.pyexammesh_faces, vertices=math5.pyexammesh_vertices,
               style="matlab")

#rotation and update of vertices
n11=0.05*np.pi
n22=0.05*np.pi

#Scale and rotate mesh
MHexam.scale_mesh(scale1)
MHexam.rotate_mesh(n11, "xy")    
MHexam.rotate_mesh(n22,"zx")

#Initialisation of edges
MHexam.unique_edges()
#Identification of edges
MHexam.silhouette_vertices("Lodox", ViewP)
#projecting alignment points
points_proj=pm.project_points_ss(MHexam.silvert1, l0=ViewP, p0=p0) 
        
#PLOT 1: Plotting mesh, silhouette points and projected points
print "Plot figure z-axis reversed (because z=0 is top of scan)"
mlab.figure(1, bgcolor=(1,1,1), size=(800,800))
mlab.clf()
mlab.triangular_mesh(MHexam.vertices[:,0], MHexam.vertices[:,1], 
                             MHexam.vertices[:,2]*(-1), MHexam.faces,
                             colormap='Greens', representation='surface',
                             opacity=0.2)

#points on the airway
mlab.points3d(MHexam.silvert1[:,0], MHexam.silvert1[:,1], MHexam.silvert1[:,2]*(-1), 
              scale_factor=1, color=(0,0,0))
#    #projected points
mlab.points3d(points_proj[:,0], points_proj[:,1], points_proj[:,2]*(-1), 
              scale_factor=1, color=(0,0,0))  

# Save figure as x3d or png
#mlab.savefig("figure1.x3d") #exporting as x3d
#mlab.savefig("figure1.png", magnification=2)
mlab.show()
