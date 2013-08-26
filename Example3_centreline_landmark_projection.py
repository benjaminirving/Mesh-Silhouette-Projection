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
proj_type="Lodox"#scipy.nd.filters.convolve

# reading in mesh and statistical model
math5=H5fileData("data_model/meshdata.h5")
scale1=1.0/math5.scale1.mean()

#Full example mesh
MHexam=mc.Mesh(faces=math5.pyexammesh_faces, vertices=math5.pyexammesh_vertices,
               style="matlab")
#PCA operations on the mesh
MH=mc.PcaMesh(faces=math5.pytempmesh_faces, vertices=math5.pytempmesh_vertices,
           style="matlab")

# P - Eigenvectors
# V - Contribution of eigenvalues
MH.input_pca(P=math5.P, V=math5.V, points=math5.pyairpoints, skel=math5.pyairskel)

###################################################################################
#3D mesh manipulation and silhouette detection

#rotation and update of vertices
n11=0.05*np.pi
n22=0.05*np.pi

#3D mesh
MH.scale_mesh(scale1)
MHexam.scale_mesh(scale1)
MH.rotate_mesh(n11, "xy")
MHexam.rotate_mesh(n11, "xy")    
MH.rotate_mesh(n22, "zx")
MHexam.rotate_mesh(n22,"zx")
#Initialisation of edges
MH.unique_edges()
#Identification of edges
MH.silhouette_vertices("Lodox", ViewP)
#projecting alignment points
points_bif3Dproj=pm.project_points_ss(MH.points, l0=ViewP, p0=p0) 
skel_3Dproj=pm.project_points_ss(MH.skel, l0=ViewP, p0=p0)
points_proj=pm.project_points_ss(MH.silvert1, l0=ViewP, p0=p0) 
        
#PLOT 1: Plotting
meshplot_module.plot_big(MH, MHexam, points_proj, points_bif3Dproj, skel_3Dproj, ViewP, proj_type, fnum=1)
#mlab.savefig("figure1.x3d") #exporting as x3d
#mlab.savefig("figure1.png", magnification=2)
mlab.show()

