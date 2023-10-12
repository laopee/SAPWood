# This code file contains major classes in SAPWood analysis
import numpy as np
import os as os
import math
from tkinter import ttk


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
   
    def EstimateK(self,new_X):
        # this function will nudge the current spring with a new location X, and find the secant stiffness K to get there
        # dont touch anything, dont update any variables, just hypothetically go to X
        new_F=self.GetNewForce(new_X)
        return (new_F-self.CuF)/(new_X-self.CuX)

    def Push(self, new_X,new_V,new_A,new_F,new_K,new_tracker):
        
        self.PaA=self.CuA
        self.PaV=self.CuV
        self.PaX=self.CuX
        self.PaF=self.CuF
        self.PaK=self.CuK   #literally pushing spring one step forward
        
        # put the New Current values in Spr status
        self.CuX=new_X
        self.CuF=new_F
        self.CuA=new_A
        self.CuV=new_V
        self.CuK=new_K  
        
        #update max-min
        self.Xmax=max(self.Xmax,self.CuX)
        self.Fmax=max(self.Fmax,self.CuF)
        self.Xmin=min(self.Xmin,self.CuX)
        self.Fmin=min(self.Fmin,self.CuF)
        self.tracker=new_tracker

        # store current values into the history
        # self.WriteCurrent()  # we dont do this here because we want to have the option of NOT record spring history to save time and memory
     
    
    
    def WriteCurrent(self):
        self.A=np.append(self.A,self.CuA)
        self.V=np.append(self.V,self.CuV)
        self.X=np.append(self.X,self.CuX)
        self.F=np.append(self.F,self.CuF)
        self.K=np.append(self.K,self.CuK)
        self.CurrentIndx+=1  # you will know how many step you are at in recording

    def ClearMemory(self):
        self.X=[]
        self.V=[]
        self.A=[]
        self.F=[]
        self.K=[]    # empty out history, but dont initialize tracker and max
        self.CurrentIndx=0
    
    def Protocal_Push(self,pro,scale):  #just do disp-control push of a spring
        
        self.__init__()
        self.init_tracker()
        
        xx=pro.value*scale   #pro is protocol object
        for ii in range(0,len(xx)):
            tempF=self.GetNewForce(xx[ii])
            tempK=self.EstimateK(xx[ii])
            tempT=self.Estimate_tracker(xx[ii])

            self.Push(xx[ii],0,0,tempF,tempK,tempT)
    
    #following are functions passed to specific class

    def SetParameter(self,inputP):
        pass  # this will be different for each spring type

    def init_tracker(self):
        pass
    
    def Estimate_tracker(self,new_X): # this function will return tracker if the spring goes to new_X
        # but it will NOT actually update or touch tracker, just hypothetically go to X
        pass
    def GetNewForce(self,new_X):   #this is the difference in different Springs
        pass                    # get force will NOT push the spring, it will just get force but no updates on variable        
    def GetK0(self):            #obtain initial K
        pass


# Assign spring based on ID
def Assign_Spr(ID:int)-> Spring:
    if ID==1:
        return Spr_Linear()
    if ID==2:
        return Spr_Bilinear()

# Spring sub class Spr_Linear  ID=1
class Spr_Linear(Spring):
    # linear spring has only 1 parameter k
    def SetParameter(self, inputP):
        self.parameter=inputP
        self.CuK=self.parameter[0]  #set initial k as k
        
    def init_tracker(self):
        self.tracker=0
        
    def GetNewForce(self,new_X):
        return new_X*self.parameter[0]
    
    def Estimate_tracker(self,new_X):
        return 0
    
    def GetK0(self):
        return self.parameter[0]

# Spring sub class Spr_Bilinear   ID=2
class Spr_Bilinear(Spring):
    # linear spring has only 3 parameter [k0 ky Dy]
    def init_tracker(self):
        self.tracker=0      # tracker is X0, the balanced location

    def GetNewForce(self,new_X):
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
        
    def GetK0(self):
        return self.parameter[0]
# Spring sub class Spr_CUREE   ID=3

# Spring sub class Spr_EPHM   ID=4

# Spring sub class Spr_MultiLinear    ID=5

# Spring sub class Spr_CompOnly    ID=6

# Spring sub class Spr_TensionOnly   ID=7

# Model_file super class
class Model_file:
    def __init__(self) -> None:
        pass        # create internal variables to store model info

    def LoadFile(self,fileLoc):
        pass        # load from file

    def SaveFile(self,fileLoc):
        pass        # save to a file

    def To_str(self)-> str:
        pass

class Model_file_SDOF(Model_file):
    
    def __init__(self) -> None:
        self.mass=0
        self.spr_type=0
        self.spr_parameter=[]

    def LoadFile(self, fileLoc):
        with open(fileLoc,'r') as file:
            lines=file.readlines()
        
        self.mass=float(lines[0].strip())
        self.spr_type=float(lines[1].strip())
        temp=lines[2].strip()
        self.spr_parameter=[float(num) for num in temp.split()]

        #SDOF file format
        #1 mass
        #2 Sprtype 1-linear  2-bilinear
        #3 SprParameters k0. Or  k0 ky Dy

    def SaveFile(self, fileLoc):
        with open(fileLoc, 'w') as file:
            # Write the first number (float) to the first line
            file.write(f"{self.mass:.6f}\n")

            # Write the second number (integer) to the second line
            file.write(f"{int(self.spr_type)}\n")

            # Write the array of numbers to the third line (space-separated)
            array_str = " ".join(str(num) for num in self.spr_parameter)
            file.write(array_str)
            
    def To_str(self) -> str:
        tempstr=""
        tempstr+=str(self.mass)+"\n"
        tempstr+=str(self.spr_type)+"\n"
        tempstr+="\t".join(str(f) for f in self.spr_parameter)

        return tempstr

# Model_Dyn super class Model
class Model_Dyn:

    def __init__(self) -> None:
        pass        # create variable spaces

    def Construct(self,Modelfile:Model_file):
        pass        # construct variables based on model

    def Initialize(self):
        pass        # renew and make the model "new"

    def Analysis_NB(self, EQ:Earthquake,ScaleFactor,DampR, timeStep, Pro_Bar, WriteSpringOK):
        pass        # earthquake analysis

    def Analysis_Push(self,Pro:Protocols,DOF_ID,ScaleFactor):
        pass        # let's try this, it is a SDOF displacement-control push

    def DofPlot(self,DOF_ID,TargetCanvas):
        pass        # plot any DOF over time

    def HystPlot(self,Spr_ID,TargetCanvas):
        pass        # plot any spring element's hystersis


# Model sub class Model_SDOF
class Model_Dyn_SDOF(Model_Dyn):
    def __init__(self):
        self.time=[]
        self.Spr=Spring()  # we don't know what spring type here, do that in construct
        self.mass=0
        self.dampR=0
        self.GlobalX=[]
        self.GlobalV=[]
        self.GlobalA=[]
        self.CurrentIndex=0

        self.CuX=0
        self.CuV=0
        self.CuA=0

    def Construct(self, Modelfile: Model_file_SDOF):
        self.mass=Modelfile.mass
        self.Spr=Assign_Spr(Modelfile.spr_type)
        self.Spr.SetParameter(Modelfile.spr_parameter)

    def Initialize(self):
        self.CurrentIndex=0
        self.GlobalX=[0]
        self.GlobalV=[0]
        self.GlobalA=[0]
        self.time=[0]
        self.Spr.ClearMemory()

        self.CuX=0
        self.CuV=0
        self.CuA=0

        self.PaX=0
        self.PaV=0
        self.PaA=0

    def Analysis_NB(self, EQ: Earthquake, ScaleFactor, DampR, timeStep, Pro_Bar:ttk.Progressbar, WriteSpringOK):
        # basic parameter
        Beta=1/6
        damping=2*math.sqrt(self.mass*self.Spr.CuK)*DampR
        # remap ground motion into timestep
        Tmax=max(EQ.t)
        tt=np.arange(0,Tmax,timeStep)
        Ft=self.mass*EQ.Ax*ScaleFactor
        Ftt=np.interp(tt,EQ.t,Ft)

        nn=len(tt)
        # set up progress bar
        Pro_Bar["maximum"]=nn
        for ii in range(len(tt)-1):
            Df=Ftt[ii+1]-Ftt[ii]
            Df_b=Df+self.CuV*(self.mass/timeStep/Beta+damping/2/Beta)+self.CuA*(self.mass/2/Beta-damping*timeStep*(1-1/4/Beta))
            K_b=self.Spr.CuK+self.mass/Beta/timeStep/timeStep+damping/2/Beta/timeStep

            D_x=Df_b/K_b
            D_v=D_x/2/Beta/timeStep-self.CuV/2/Beta+self.CuA*timeStep*(1-1/4/Beta)

            #update model X and V
            self.PaX=self.CuX
            self.PaV=self.CuV
            self.CuX=self.PaX+D_x
            self.CuV=self.PaV+D_v


            #for all springs in the model calculate new properties
            # potentially you can implement sub-step here
            temp_newTrack=self.Spr.Estimate_tracker(self.CuX)   #get new tracker
            temp_newF=self.Spr.GetNewForce(self.CuX)  #get new spring force
            temp_newK=self.Spr.EstimateK(self.CuX)      #get new spring K
            
            #with spring new forces, update model A
            self.PaA=self.CuA
            self.CuA=(Ftt[ii+1]-damping*self.CuV-temp_newF)/self.mass

            #update Springs
            self.Spr.Push(self.CuX,self.CuV,self.CuA,temp_newF,temp_newK,temp_newTrack)
            if WriteSpringOK:  #write Spring data if desired
                self.Spr.WriteCurrent()
            
            #write in model vectors
            self.GlobalX=np.append(self.GlobalX,self.CuX)
            self.GlobalV=np.append(self.GlobalV,self.CuV)
            self.GlobalA=np.append(self.GlobalA,self.CuA)
            self.time=np.append(self.time,tt[ii+1])

            #progress bar
            Pro_Bar["value"]=ii
            Pro_Bar.update()
            #might need to get root or pbar to update




     


# Model sub class Model_ShearBld

# Model sub class Model_Biaxial

# Model sub class Model_6DOF
