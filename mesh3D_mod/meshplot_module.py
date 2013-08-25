"""
Using this file to test some plotting functions. 
"""

import numpy as np
import matplotlib.pyplot as plt
import mayavi.mlab as mlab


def test_plot3d():
    """Generates a pretty set of lines."""
    n_mer, n_long = 6, 11
    pi = np.pi
    dphi = pi/1000.0
    phi = np.arange(0.0, 2*pi + 0.5*dphi, dphi)
    mu = phi*n_mer
    x = np.cos(mu)*(1+np.cos(n_long*mu/n_mer)*0.5)*100
    y = np.sin(mu)*(1+np.cos(n_long*mu/n_mer)*0.5)*100
    z = np.sin(n_long*mu/n_mer)*0.5*100

    l = mlab.points3d(x, y, z, scale_factor=5)
    return l
    
def plot_big(MH, MHexam, points_proj, points_projal, skel_projal, ViewP, proj_type, plot_type="normal", fnum=1):
  
    # Plotting
    print "Plot figure z-axis reversed (because z=0 is top of scan)"
    mlab.figure(fnum, bgcolor=(1,1,1), size=(800,800))
    mlab.clf()
    #airway plot    
    mlab.triangular_mesh(MH.vertices[:,0], MH.vertices[:,1], MH.vertices[:,2]*(-1), 
                         MH.faces, colormap='Blues', representation='wireframe', 
                         line_width=0.5)
    try:
        mlab.triangular_mesh(MHexam.vertices[:,0], MHexam.vertices[:,1], 
                             MHexam.vertices[:,2]*(-1), MHexam.faces,
                             colormap='Greens', representation='surface',
                             opacity=0.2)
    except:
        print "Example mesh not plotted"
    #airway axes
    #mlab.axes('off') #extent=([-200, 200, -200,200, 200,-200])) 
    #points on the airway
    mlab.points3d(MH.silvert1[:,0], MH.silvert1[:,1], MH.silvert1[:,2]*(-1), 
                  scale_factor=1, color=(0,0,0))
#    #projected points
    mlab.points3d(points_proj[:,0], points_proj[:,1], points_proj[:,2]*(-1), 
                  scale_factor=1, color=(0,0,0))  
    
    if plot_type is "normal":        
#    #alignment points
        mlab.points3d(MH.points[:,0], MH.points[:,1], MH.points[:,2]*(-1), 
                      scale_factor=2, color=(0,1,0))
        
        #skeleton
        mlab.points3d(MH.skel[:,0], MH.skel[:,1], MH.skel[:,2]*(-1), 
                      scale_factor=0.5, color=(0.5,0.5,1))
    #    #alignment points
        mlab.points3d(points_projal[:,0], points_projal[:,1], points_projal[:,2]*(-1), 
                      scale_factor=2, color=(0,1,0))
        
    try:
        mlab.points3d(skel_projal[:,0], skel_projal[:,1], skel_projal[:,2]*(-1), 
                      scale_factor=0.5, color=(0.5,0.5,1))
    except:
        print "skeleton not plotted"
    #view direction
    #if proj_type=="Lodox":
    #    ViewL=np.array([ViewP, ViewP])
    #    ViewL[0,2]=30; ViewL[1,2]=-30        
    #    
    #    mlab.plot3d(ViewL[:,0]/10, ViewL[:,1], ViewL[:,2]*(-1), 
    #                color=(0,0,0), 
    #                tube_radius=1, opacity=1)
    
    #elif proj_type=="normal":
    #    mlab.points3d(ViewP[0]/10, ViewP[1], ViewP[2]*(-1), 
    #                  scale_factor=4, color=(1,0,0))
    
    mlab.view(elevation=80,azimuth=20)
    #fname="%s/figure_save2/tmp%04i.png" %(os.getcwd(), cc1);
    #mlab.savefig(fname)

def plot_im(Im, Pnt):  
    
    plt.figure(3)
    plt.imshow(Im.I, cmap=plt.cm.gray)
    #plt.axis([0, Im.I.shape[1], 0, Im.I.shape[0]])
    #plot in reverse plot is for (x,y) not (rows, columns)
    plt.plot(Im.yy, Im.xx, 'ro')
    plt.plot(Im.yy[3], Im.xx[3], 'r*')
    #plot in reverse plot is for (x,y) not (rows, columns)
    plt.plot(Pnt.match3d_align[:,1], Pnt.match3d_align[:,0], 'go')
    plt.plot(Pnt.match3d_align[3,1], Pnt.match3d_align[3,0], 'g*')
    plt.plot(Pnt.point2d[:,1], Pnt.point2d[:,0], 'g.')