# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 15:38:21 2012

@author: b1e5n98jamin
"""
import h5py

class H5fileData:
    
    #want everything to happen in initialisation
    def __init__(self, h5file):
        self.h5file=h5file
        self.__load()
    
    # loading files (but used in init)
    def __load(self):
        f=h5py.File(self.h5file, 'r')
        #vertices=f['/vertices'].value
        
        for ii in list(f):
#            print ii
            setattr(self, ii, f[ii].value.transpose())
            