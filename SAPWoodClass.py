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
        ncol=temp.shape[1]
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
class Protocols: 
    # All protocol obj has a maximum value of 1, so all need to be scaled before using
    

    def __init__(self):
        self.step_size=0.01
        self.value=np.arange(0,1,self.step_size) # by default this will be a monotonic protocol from 0 to 1 with given step size
        self.max=1
        
    def changeStepSize(self,newsize):
        self.step_size=newsize

    def scale(self,targetMax):
        factor=targetMax/self.max
        self.value=np.multiply(self.value,factor)
        self.max=targetMax

    def cyclic_linear(self,N_cycle): # Linear increasing amplitude with step size
        keypt=np.arange(0,1/N_cycle,1)
        ii=0
        s=-1
        e=1
        self.value=[]
        while ii<len(keypt):
            temp=np.arange(s*keypt[ii],self.step_size,e*keypt[ii+1])
            s=-s
            e=-e
            ii+=1
            self.value=np.append(self.value,temp)

# General loading class Load

# Spring super class Spr
class Spring:
    def __init__(self):
        self.parameter=[] # spring parameters
        self.tracker=[]  # special tracker info for nonlinear springs
        
        self.X=[]
        self.V=[]
        self.A=[]
        self.F=[]
        self.K=[]    # these are records of spring history
        self.Xmax=0
        self.Fmax=0
        self.Xmin=0
        self.Fmin=0  #max and min
        self.CuX=0
        self.CuV=0
        self.CuA=0
        self.CuK=0
        self.CuF=0  # current spring values

        self.PaX=0
        self.PaV=0
        self.PaA=0
        self.PaK=0
        self.PaF=0  # past spring values (last step)

        self.CurrentIndx=0  # where spring is currently

    def SetParameter(self,inputP):
        self.parameter=inputP  # this will be different for each spring type

    def init_tracker(self):
        pass
    
    def Estimate_tracker(self,new_X): # this function will return tracker if the spring goes to new_X
        # but it will NOT actually update or touch tracker, just hypothetically go to X
        pass

    def EstimateK(self,new_X):
        # this function will nudge the current spring with a new location X, and find the secant stiffness K to get there
        # dont touch anything, dont update any variables, just hypothetically go to X
        new_F=self.GetForce(new_X)
        return (new_F-self.CuF)/(new_X-self.CuX)

    def Push(self, new_X,new_V,new_A,new_F,new_K,new_tracker):
        
        self.PaA=self.CuA
        self.PaV=self.CuV
        self.PaX=self.CuX
        self.PaF=self.CuF
        self.PaK=self.CuK   #literally pushing spring one step forward
        
        # calculate the New Current values
        self.CuX=new_X
        self.CuF=new_F
        self.CuA=new_A
        self.CuV=new_V
        self.CuK=new_K  #this functionwill find force
        
        #update max-min
        self.Xmax=max(self.Xmax,self.CuX)
        self.Fmax=max(self.Fmax,self.CuF)
        self.Xmin=min(self.Xmin,self.CuX)
        self.Fmin=min(self.Fmin,self.CuF)
        self.tracker=new_tracker

        # store current values into the history
        self.WriteCurrent()
        self.CurrentIndx+=1
    
    def GetForce(self,new_X):   #this is the difference in different Springs
        pass                     # get force will NOT push the spring, it will just get force but no updates on variable        

    def WriteCurrent(self):
        self.A=np.append(self.A,self.CuA)
        self.V=np.append(self.V,self.CuV)
        self.X=np.append(self.X,self.CuX)
        self.F=np.append(self.F,self.CuF)
        self.K=np.append(self.K,self.CuK)

    def ClearMemory(self):
        self.X=[]
        self.V=[]
        self.A=[]
        self.F=[]
        self.K=[]    # empty out history, but dont initialize tracker and max
        self.CurrentIndx=0
    
    def Protocal_Push(self,pro,scale):
        
        self.__init__()
        self.init_tracker()
        
        xx=pro.value*scale   #pro is protocol object
        for ii in range(0,len(xx)):
            tempF=self.GetForce(xx[ii])
            tempK=self.EstimateK(xx[ii])
            tempT=self.Estimate_tracker(xx[ii])

            self.Push(xx[ii],0,0,tempF,tempK,tempT)


# Spring sub class Spr_Linear
class Spr_Linear(Spring):
    # linear spring has only 1 parameter k
    def init_tracker(self):
        self.tracker=0
        
    def GetForce(self,new_X):
        return new_X*self.parameter[0]
    
    def Estimate_tracker(self,new_X):
        return 0

# Spring sub class Spr_Bilinear
class Spr_Bilinear(Spring):
    # linear spring has only 3 parameter [k0 ky Dy]
    def init_tracker(self):
        self.tracker=0      # tracker is X0, the balanced location

    def GetForce(self,new_X):
        K0=self.parameter[0]
        Ky=self.parameter[1]
        Dy=self.parameter[2]
        X0=self.tracker

        if abs(self.CuX-X0)<=Dy:  # if currently we are in K0 region
            if abs(new_X-X0)<=Dy:
                return X0*Ky+(new_X-X0)*K0
            elif new_X-X0>0:
                return X0*Ky+Dy*K0+(new_X-X0-Dy)*Ky
            elif new_X-X0<0:
                return X0*Ky-Dy*K0-(X0-new_X-Dy)*Ky
            else:
                print("sth wrong with Bilinear spring")

        if self.CuX-X0>Dy:   # if we are in positive yielding region
            if new_X-self.CuX>0:# and going more positive
                return X0*Ky+Dy*K0+(new_X-X0-Dy)*Ky
            else:                   # and going negative
                # could update tracker here
                if new_X-self.CuX>-Dy*2:    # go negative a little
                    return self.CuF+K0*(new_X-self.CuX)
                else:                       # go negative a lot
                    return self.CuF-K0*2*Dy+Ky*(new_X-self.CuX+2*Dy)
        
        if self.CuX-X0<-Dy:   # if we are in negative yielding region
            if new_X-self.CuX<0:   # and going more negative
                return X0*Ky-Dy*K0-(X0-new_X-Dy)*Ky
            else:                   # and going positive
                # could update tracker here
                if new_X-self.CuX<Dy*2:  # positive a little
                    return self.CuF+K0*(new_X-self.CuX)
                else:                   #  positive a lot
                    return self.CuF+K0*2*Dy+Ky*(new_X-self.CuX-2*Dy)
    
    def Estimate_tracker(self, new_X):
        K0=self.parameter[0]
        Ky=self.parameter[1]
        Dy=self.parameter[2]
        X0=self.tracker

        if abs(self.CuX-X0)<=Dy:  # if currently we are in K0 region
            return X0               # no change in X0, so return the same X0

        if self.CuX-X0>Dy:   # if we are in positive yielding region
            if new_X-self.CuX>0:# and going more positive
                return X0
            else:                   # and going negative
                # could update tracker here
                return self.Cux-Dy
        
        if self.CuX-X0<-Dy:   # if we are in negative yielding region
            if new_X-self.CuX<0:   # and going more negative
                return X0
            else:                   # and going positive
                return self.CuX+Dy
        
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
