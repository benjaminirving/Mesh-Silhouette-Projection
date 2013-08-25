""" Point manipulation module
Functions on this module are used for projection, transformation and manipulation of points in
3D and 2D space. 
"""

from numpy import dot, array, sqrt
import numpy as np
import matplotlib.pyplot as plt


def project_points(points1, p0=[-400.0, 0.0, 0.0], n=[1.0,0,0], l0=[100.0, 0.0, 0.0]):
    """
    Project silhouette points onto a plane from a point source
    
    p0 is point on the plane
    n is the plane normal
    l0 is point from view
    points are the points that are projected
    """
    
    p0=array(p0); n=array(n); l0=array(l0);
    
    l1=points1-l0
    
    l1_norm=((points1**2).sum(axis=1))**0.5
    l1_norm=l1_norm.reshape(l1_norm.shape[0], 1)
    l1=l1/l1_norm
    
    #solving for d for each vertex (d=((p0-l0) dot n) / (l dot n) )
    d1=dot((p0-l0), n)/((l1*n).sum(axis=1))
    d1=d1.reshape(d1.shape[0],1)
    
    #projected points
    points_proj=d1*l1 + l0.reshape(1,3)
    
    return points_proj


def project_points_ss(points1, p0=[-400.0, 0.0, 0.0], n=[1.0,0,0], l0=[100.0, 0.0, 0.0]):
    """
    Project silhouette points onto a plane using a linear slit scanning projection
    
    p0 is point on the plane
    n is the plane normal
    l0 is point from view
    points are the points that are projected
    """
    
    p0=array(p0); n=array(n); l0=array(l0);
    
    #creating l0 with the same z values as the points about to be projected
    l0proc=np.tile(l0, (points1.shape[0], 1))
    l0proc[:,2]=points1[:,2]    
    
    #vector from viewpoint to points on mesh
    l1=points1-l0proc
    #normalise
    l1_norm=((points1**2).sum(axis=1))**0.5
    l1_norm=l1_norm.reshape(l1_norm.shape[0], 1)
    l1=l1/l1_norm
    
    #solving for d for each vertex (d=((p0-l0) dot n) / (l dot n) )
    d1=dot((p0-l0), n)/((l1*n).sum(axis=1))
    d1=d1.reshape(d1.shape[0],1)
    
    #projected points
    points_proj=d1*l1 + l0proc
    
    return points_proj
    
def similarity_transform2D(p1, p2):
    """ 
    translation, scaling and rotation of p1 onto p2
    (see Stegman2002bis - A brief introduction to statistical shape analysis
    as well as "kheng000rr - rigid registration - rigid - similarity - affine.pdf")            
    """
    p1orig=p1
    p2orig=p2
    
    #Translation (translate both sets of vectors to their centroid)
    m1=p1.mean(axis=0)
    m2=p2.mean(axis=0)
    p1=p1-m1
    p2=p2-m2
    #plt.figure(1)
    #plt.plot(p1[:,0], p1[:,1],'ro', p2[:,0], p2[:,1], 'go')
    
    #Scale p1 to p2 (Frobenius norm)
    p1s=sqrt(sum(sum(p1**2)))
    p2s=sqrt(sum(sum(p2**2)))
    s=p2s/p1s;
    p1=p1*s;
    #plt.figure(2)
    #plt.plot(p1[:,0], p1[:,1],'ro', p2[:,0], p2[:,1], 'go')
    
    #Rotation (SVD of correlation to find optimal rotatin of p1)
    corr1=np.dot(p1.transpose(), p2)
    u, ss, v=np.linalg.svd(corr1)
    R=np.dot(v, u.transpose())
    p1=np.dot(R, p1.transpose())
    p1=p1.transpose()    
    #plt.figure(3)
    #plt.plot(p1[:,0], p1[:,1],'ro', p2[:,0], p2[:,1], 'go')
    
    #Calculate transform T
    m1prime=s * np.dot(R, m1.transpose())
    m1prime=m1prime.transpose()
    T=m2-m1prime
    #p1orig2=s*(np.dot(R, p1orig.transpose())).transpose() + T
    #plt.figure(4)
    #plt.plot(p1orig2[:,0], p1orig2[:,1],'ro', p2orig[:,0], p2orig[:,1], 'go')
    #plt.show()
    
    return s, R, T

        
