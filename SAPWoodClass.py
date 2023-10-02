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
        temp=np.genfromtxt(filename,delimiter=' ', names=['c1','c2','c3','c4'])
        self.t=temp['c1']
        self.Ax=temp['c2']
        self.Ay=temp['c3']
        self.Az=temp['c4']

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
