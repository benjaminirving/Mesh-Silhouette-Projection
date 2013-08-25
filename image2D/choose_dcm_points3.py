#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 17:42:54 2012

This version saves image as a temporary .jpg file. This bypasses the use
of qimage2ndarry that seems to have issue

run the file from the command line:
python choose_dcm_points3.py dir <directory>
python choose_dcm_points3.py <type> <dicom_name> 

dicom_name: dicom file name including '.dcm'
type: 'align', 'edge' or put nothing. Appends word to saved file
directory: directory where file is found. Only neccessary for first file

choose directory first:
python choose_dcm_points3.py dir /home/b1e5n98jamin/CADdata/LodoxTBanon/dicoms/

select alignment points:
python choose_dcm_points3.py align 'B03 -APLo' 

select edge points:
python choose_dcm_points3.py edge 'B03 -APLo'

sometimes it is neccessary to have multiple measurements
(if a 3rd value it appended then it will be added to the output file)

python choose_dcm_points3.py edge 'B03 -APLo' 2
    
@author: benjamin
"""



import sys

from PyQt4 import QtCore, QtGui, Qt
import numpy as np

class MainWidget(QtGui.QWidget):
    
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.scene = QtGui.QGraphicsScene()
        self.view = QtGui.QGraphicsView(self.scene)
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.view)
        self.setLayout(layout)
        self.pixmap_item = QtGui.QGraphicsPixmapItem(QtGui.QPixmap('tmp129495324.jpg'), None, self.scene)
        self.pixmap_item.mousePressEvent = self.pixelSelect
        self.xx=[]
        self.yy=[]


    def pixelSelect(self, event):
        x=event.pos().y();
        y=event.pos().x();
        print "Click!", x, y
        self.xx.append(x)   
        self.yy.append(y)
        
        pen = QtGui.QPen(QtCore.Qt.red)
        if len(self.xx)>1:
            self.line1=Qt.QLineF(self.yy[-2], self.xx[-2], self.yy[-1], self.xx[-1])
            self.scene.addLine(self.line1, pen)
            
            #self.scene.addPolygon(QtGui.QPolygonF(self.line1), pen)
        #for ii in range(len(self.xx)):
        self.scene.addEllipse(y-5, x-5, 10, 10, pen)
        
    def return_vals(self):
        coords1=np.zeros((len(self.xx),2))
        coords1[:,0]=self.xx
        coords1[:,1]=self.yy
        return coords1
        


if __name__ == "__main__":
    
    #import dicom
    import imagefilt2D_module as ft
    #import qimage2ndarray as qnd
    import pickle
    import scipy as sp
    
    #Looking for command line inputs
  
        
    # if the file directory is included 
    #(save directory and exit)
    
    if sys.argv[1]=="help" or sys.argv[1]=="--help" or sys.argv[1]=="-h":
        print __doc__
        sys.exit()
    
    if sys.argv[1]=="dir":
        try:
            print sys.argv[2]
            filedir1=sys.argv[2]
            
            # save file directory
            dinp=open("tmp23423423.dat", 'w')
            pickle.dump(filedir1, dinp)
            dinp.close()
            print "Done"
        except:            
            sys.exit("Directory has NOT been stored")
        sys.exit("Directory has been stored")
    # if align or edge then load image 

    elif sys.argv[1]=="align":
        print "Alignment points"
        print "Trachea, BIF, LMB, RMB"
        name_ad="_align"
    elif sys.argv[1]=="edge":
        print "Points around the edge"
        name_ad="_edge"
    elif sys.argv[1]=="test":
        name_ad="_test"
    else:
        sys.exit("No such command")
    
    try:
        dinp=open("tmp23423423.dat", 'r')
        filedir1=pickle.load(dinp)
        dinp.close()
    except:
        sys.exit("No file stored!!!")
        
    #extra appended file name
    try:
        extran=sys.argv[3]
    except:
        extran=""
        
    fname=sys.argv[2]
    print fname[0:-4]
    
    Im=ft.XrayImage(filedir1, fname)
    Im.crop_border()
    
    Im.doub_norm()
    Im.I=ft.local_norm(Im.I,20)
    
    sp.misc.imsave('tmp129495324.jpg', Im.I)
#    
#    #Exporting to qimage
    #qimage1 = qnd.gray2qimage(Im.I, (Im.I.min(), Im.I.max()))
    
    app = QtGui.QApplication(sys.argv)
    widget = MainWidget()
    widget.resize(900, 700)
    widget.show()
    app.exec_()
    
    xx =widget.return_vals()
    print xx
    
    if name_ad=="_align" and xx.shape[0] != 4:
        sys.exit("Error! wrong number of points for align")

    if name_ad=="_edge" and xx.shape[0] <= 4:
        sys.exit("Error! wrong number of points for edge")
    
    save_y=raw_input("Save data? y/n: ")
    
    if save_y=="y" or save_y=="Y":
        print "Saving..."
        output=open(filedir1 + fname + name_ad + extran + '.pkl', 'wb')
        pickle.dump(xx, output)
        output.close()
        print "Done"


