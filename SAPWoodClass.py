# This code file contains major classes in SAPWood analysis
import numpy as np
import os as os


# Earthquake class EQ
class Earthquake:

    def _init_(self):
        self.name='Null'
        self.Ax=np.array()
        self.Ay=np.array()
        self.Az=np.array()
        self.t=np.array()
        self.dir=3 #default to have all 3 directions
        self.Unit='G' #default unit is in G 


    def LoadEQ(self,filename):  #this loads earthquake from a text file
        temp=np.loadtxt(filename)
        ncol=temp.shap[1]
        if ncol==2:
            self.t=temp[:,0]
            self.Ax=temp[:,1]
            self.Ay=np.copy(self.Ax)  #fill other directions with 0
            self.Ay.fill(0)
            self.Az=np.copy(self.Ay)

        elif ncol==3:
            self.t=temp[:,0]
            self.Ax=temp[:,1]
            self.Ay=temp[:,2]
            self.Az=np.copy(self.Ay)
            self.Az.fill(0)
        elif ncol==4:
            self.t=temp[:,0]
            self.Ax=temp[:,1]
            self.Ay=temp[:,2]
            self.Az=temp[:,3]
        else:
            print("something wrong with EQ file, pls check")

    def GetDt(self):
        return self.t[1]-self.t[0]



# Push-over class Push

# Displacement protocol class Disp

# General loading class Load

# Spring super class Spr

# Spring sub class Spr_Linear

# Spring sub class Spr_Bilinear

# Spring sub class Spr_CUREE

# Spring sub class Spr_EPHM

# Spring sub class Spr_MultiLinear

# Spring sub class Spr_CompOnly

# Spring sub class Spr_TensionOnly

# Model super class Model

# Model sub class Model_SDOF

# Model sub class Model_ShearBld

# Model sub class Model_Biaxial

# Model sub class Model_6DOF
