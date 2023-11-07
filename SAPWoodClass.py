# This code file contains major classes in SAPWood analysis
import numpy as np
import os as os
import math as math
from tkinter import ttk

# Utility functions used in classes
def plot_on_canvas(X, Y, canvas):  
    """
    Plots paired vectors X and Y on a canvas.

    Args:
        X (list): List of X-coordinates.
        Y (list): List of Y-coordinates.
        canvas (tk.Canvas): The canvas object to draw on.
    """
    # Determine the canvas dimensions
    canvas_width = canvas.winfo_reqwidth()
    canvas_height = canvas.winfo_reqheight()

    # Find the minimum and maximum values for X and Y
    x_min, x_max = min(X), max(X)
    y_min, y_max = min(Y), max(Y)

    print('limits:',x_min,x_max,y_min,y_max)
    print('length:',len(X),len(Y))
    #print(Y)

    # Scale the data to fit within the canvas
    x_scale = canvas_width / (x_max - x_min)
    y_scale = canvas_height / (y_max - y_min)

    # Create a blue line connecting the data points
    for i in range(len(X) - 1):
        x1 = (X[i] - x_min) * x_scale
        y1 = canvas_height - (Y[i] - y_min) * y_scale
        x2 = (X[i + 1] - x_min) * x_scale
        y2 = canvas_height - (Y[i + 1] - y_min) * y_scale
        canvas.create_line(x1, y1, x2, y2, fill="blue", width=2)

    # Mark the X and Y axes in red
    canvas.create_line(0, canvas_height - (0 - y_min) * y_scale, canvas_width, canvas_height - (0 - y_min) * y_scale, fill="red", width=2)
    canvas.create_line((0 - x_min) * x_scale, 0, (0 - x_min) * x_scale, canvas_height, fill="red", width=2)

# Assign spring based on ID   
def Assign_Spr(ID:int):
    if ID==1:
        return Spr_Linear()
    if ID==2:
        return Spr_Bilinear()
    if ID==3:
        return Spr_CUREE()
    if ID==4:
        return Spr_EPHM()
    if ID==5:
        return Spr_Multilinear()
    if ID==6:
        return Spr_CompOnly()
    if ID==7:
        return Spr_TensionOnly()
    

def Line_interXY(k1,x1,y1,k2,x2,y2):
    if abs(k1-k2)<1e-10:
        print('parallel line dont intersect')
        return x1,y1
    else:
        x=((k1 * x1 - k2 * x2) - (y1 - y2)) / (k1 - k2)
        y=k1*(x-x1)+y1
        return x,y
    
# end of all Utility functions in this class

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
        self.value0=np.arange(0,1,self.step_size) # by default this will be a monotonic protocol from 0 to 1 with given step size
        self.max=1
        self.value=self.value0 #value0 is the original data, value is what is refined and used
        
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
    
    def LoadPro(self,filename):
        temp=np.loadtxt(filename)
        self.value0=temp
        self.max=np.max(abs(temp))  
        self.value=self.value0 #value0 is the original data, value is what is refined and used      
        #  print("something wrong with Protocol file, pls check")

    def Refine(self,N_step): # refine current protocol values with N substep
        # Calculate the step size for interpolation
        #step = 1 / (N_step+ 1)
        X=self.value0

        # Initialize an empty array for interpolated points
        X_combined = []

        # Interpolate between each pair of consecutive elements
        for i in range(len(X) - 1):
            x_start, x_end = X[i], X[i + 1]
            X_interp = np.linspace(x_start, x_end, N_step + 2,endpoint=False)
            # Combine original X values with interpolated points
            X_combined = np.append(X_combined,X_interp)

        self.value=X_combined

# General loading class Load

# Spring super class Spr
class Spring:
    def __init__(self):
        self.parameter=[] # spring parameters
        self.tracker=[]  # special tracker info for nonlinear springs
        self.type=0
        
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
        if abs(new_X-self.CuX)<1e-10:
            print("divide by 0, CuK unchanged")
            return self.CuK
        else:
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
        
        self.ClearMemory()
        self.init_tracker()
        
        xx=pro.value*scale   #pro is protocol object
        for ii in range(0,len(xx)):
            tempF=self.GetNewForce(xx[ii])
            tempK=self.EstimateK(xx[ii])
            tempT=self.Estimate_tracker(xx[ii])

            self.Push(xx[ii],0,0,tempF,tempK,tempT)
            self.WriteCurrent()
    
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

    def HysPlot(self,canvas):
        X=self.X
        Y=self.F

        plot_on_canvas(X,Y,canvas)

        # Scale the data points to fit within the canvas
        #x_min, x_max = min(X), max(X)
        #y_min, y_max = min(Y), max(Y)

        #for i in range(len(X) - 1):
        #    # Map X and Y to canvas coordinates for each pair of adjacent points
        #    x1 = (X[i] - x_min) * canvas.winfo_width() / (x_max - x_min)
        #    y1 = canvas.winfo_height() - (Y[i] - y_min) * canvas.winfo_height() / (y_max - y_min)
        #    x2 = (X[i + 1] - x_min) * canvas.winfo_width() / (x_max - x_min)
        #    y2 = canvas.winfo_height() - (Y[i + 1] - y_min) * canvas.winfo_height() / (y_max - y_min)

        #    # Draw a line segment connecting adjacent data points
        #    canvas.create_line(x1, y1, x2, y2, fill="blue")
        
        #want to add axis, maybe later. canvas.create_line(0,y_min,0,y_max,fill="red")



# Spring sub class Spr_Linear  ID=1
class Spr_Linear(Spring):
    # linear spring has only 1 parameter k
    def SetParameter(self, inputP):
        self.parameter=inputP
        self.CuK=self.parameter  #set initial k as k
        self.type=1
        
    def init_tracker(self):
        self.tracker=0
        
    def GetNewForce(self,new_X):
        #print(self.parameter[0])
        #print('inside linear get force')
        #print(self.CuK)
        return new_X*self.parameter
    
    def Estimate_tracker(self,new_X):
        return 0
    
    def GetK0(self):
        return self.parameter

# Spring sub class Spr_Bilinear   ID=2
class Spr_Bilinear(Spring):
    # linear spring has only 3 parameter [k0 ky Dy]
    def SetParameter(self, inputP):
        self.parameter=inputP
        self.CuK=self.parameter[0]  #set initial k as k0
        self.tracker=0
        self.type=2

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
                return self.CuX-Dy
        
        if self.CuX-X0<-Dy:   # if we are in negative yielding region
            if new_X-self.CuX<0:   # and going more negative
                return X0
            else:                   # and going positive
                return self.CuX+Dy
        
    def GetK0(self):
        return self.parameter[0]
# Spring sub class Spr_CUREE   ID=3
class Spr_CUREE(Spring):
    # CUREE model has 10 parameters
    # k0, xu, F0, F1, r1, r2, r3, r4, alpha, beta
    def SetParameter(self, inputP):
        self.parameter=inputP
        self.CuK=self.parameter[0]  #set initial k as k0
        self.tracker=0
        self.type=3
        self.k0=self.parameter[0]  #hard to keep writing parameter, let's do names
        self.xu=self.parameter[1]
        self.F0=self.parameter[2]
        self.F1=self.parameter[3]
        self.r1=self.parameter[4]
        self.r2=self.parameter[5]
        self.r3=self.parameter[6]
        self.r4=self.parameter[7]
        self.alpha=self.parameter[8]
        self.beta=self.parameter[9]

        self.tracker=[0, 0, 0, 0, 0, 0, 0]
        # tracker shouldbe maxX unloadX_positive minX unloadX_negative unloadX unloadF PathID
        # PathID is a number from 0-5, check the code document for meaning
        
    def init_tracker(self):
        self.tracker=[0, 0, 0, 0, 0, 0, 0]
        # tracker shouldbe maxX unloadX_positive minX unloadX_negative unloadX unloadF PathID

    def __Backbone(self,xx):
        # this is a subroutine of hysteresis model
        # find the F on the envelope curve given xx
        # xx can be + or -
        temp=0.1
        Fu = (self.F0 + self.k0 * self.r1 * self.xu) * (1 - math.exp(-self.k0 * self.xu / self.F0))
        x_end=self.xu-Fu/self.k0/self.r2
        abxx=abs(xx)

        if abxx <= self.xu:
            #print(-self.k0 * xx / self.F0)
            temp=(self.F0 + self.k0 * self.r1 * abxx) * (1 - math.exp(-self.k0 * abxx / self.F0))

        if abxx > self.xu and abxx<x_end:
            temp=Fu + (abxx - self.xu) * self.k0 * self.r2
        
        if abxx>=x_end:
            temp=0.1
        
        return np.sign(xx)*temp
    
    def __Backbone_unload_P(self,xx,unload_X):
        # this subroutine will find F on unloading path given xx
        # this only works for positive
        # you should not get this from tracker as it was not updatedunload_X=self.tracker[1]
        unload_F=self.__Backbone(unload_X)
        return unload_F+(xx-unload_X)*self.k0*self.r3
    
    def __Backbone_unload_N(self,xx,unload_X):
        # this subroutine will find F on unloading path given xx
        # this only works for negative
        
        unload_F=self.__Backbone(unload_X)
        return unload_F+(xx-unload_X)*self.k0*self.r3
    
    def __Pinch_lower(self,xx):
        # this return F on lower piching line and Kp line
        # calculated the target point on the backbone

        tar_X=self.beta*self.tracker[3]
        tar_F=self.__Backbone(tar_X)
        #print(tar_X)
        if abs(tar_X)>=1e-8:            
            kp=self.k0*(self.F0/self.k0/abs(tar_X))**self.alpha
        else:
            kp=self.k0**self.alpha
        # find the intersection of the r4 and kp line
        inX,inY=Line_interXY(kp,tar_X,tar_F,self.k0*self.r4,0,-self.F1)

        # print('intX=',inX,kp,tar_X,tar_F)
        if xx<inX and xx<-3*self.F1/self.k0:
            return tar_F+(xx-tar_X)*kp
        else:
            return xx*self.k0*self.r4-self.F1
    
    def __Pinch_upper(self,xx):
        # this return F on upper pinching line and Kp line
        # calculated the target point on the backbone
        tar_X=self.beta*self.tracker[1]
        tar_F=self.__Backbone(tar_X)
        if abs(tar_X)>=1e-8:            
            kp=self.k0*(self.F0/self.k0/abs(tar_X))**self.alpha
        else:
            kp=self.k0**self.alpha
        #kp=self.k0*(self.F0/self.k0/tar_X)**self.alpha
        # find the intersection of the r4 and kp line
        inX,inY=Line_interXY(kp,tar_X,tar_F,self.k0*self.r4,0,self.F1)

        if xx>inX and xx>3*self.F1/self.k0:
            return tar_F+(xx-tar_X)*kp
        else:
            return xx*self.k0*self.r4+self.F1
    def __General_unload(self,xx,unloadX,unloadF):
        # this return F for all unloading r3 line NOT from backbone
        # get current unload x and f and find force

        return unloadF+self.r3*self.k0*(xx-unloadX)
    

    def GetNewForce(self, new_X):
        maxX=self.tracker[0]
        unloadX_P=self.tracker[1]   #unloading form P backbone
        minX=self.tracker[2]
        unloadX_N=self.tracker[3]   #unloading from N backbone
        unloadX=self.tracker[4]     #normal unloading X
        unloadF=self.tracker[5]     #normal unloading F
        PathID=self.tracker[6]      #current location
        DX=new_X-self.CuX

        Res=0

        New_PathID=PathID

        #print('Path=',PathID,'X=',self.CuX,'Y=',self.CuF)
        #print('track=',self.tracker)

        # on positive back bone  0
        if PathID==0 and DX>=0:
            maxX=new_X
            Res=self.__Backbone(new_X)
        
        if PathID==0 and DX<0:
            unloadX_P=self.CuX
            New_PathID=4
            Res=self.__Backbone_unload_P(new_X,unloadX_P)
        
        # on negative back bone   1
        if PathID==1 and DX<=0:
            minX=new_X
            Res=self.__Backbone(new_X)
        
        if PathID==1 and DX>0:
            unloadX_N=self.CuX
            New_PathID=4
            Res=self.__Backbone_unload_N(new_X,unloadX_N)
        
        # on unloading path from backbone 4
        if PathID==4:
            if self.CuX>0:
                # find intersection of unloading from positive to pinching line
                inX,inY=Line_interXY(self.r3*self.k0,self.CuX,self.CuF,self.r4*self.k0,0,-self.F1)

                if new_X>unloadX_P:
                    maxX=new_X
                    New_PathID=0
                    Res=self.__Backbone(new_X)
                if new_X<unloadX_P and new_X>inX:
                    Res=self.__Backbone_unload_P(new_X,unloadX_P)
                if new_X<inX:
                    New_PathID=3
                    Res=self.__Pinch_lower(new_X)
            
            if self.CuX<0:
                # find intesection of unloading from negative to pinching line
                inX,inY=Line_interXY(self.r3*self.k0,self.CuX,self.CuF,self.r4*self.k0,0,self.F1)
                
                if new_X<unloadX_N:
                    minX=new_X
                    New_PathID=1
                    Res=self.__Backbone(new_X)
                if new_X>unloadX_N and new_X<inX:
                    Res=self.__Backbone_unload_N(new_X,unloadX_N)
                if new_X>inX:
                    New_PathID=2
                    Res=self.__Pinch_upper(new_X)

        # on upper r4 path 2
        if PathID==2:
            tar_X=self.beta*self.tracker[1]
            if tar_X<3*self.F1/self.k0:
                tar_X=3*self.F1/self.k0
            
            if DX>0 and new_X>=tar_X:
                maxX=new_X
                New_PathID=0
                Res=self.__Backbone(new_X)
            if DX>0 and new_X<tar_X:
                Res=self.__Pinch_upper(new_X)
            if DX<0:
                New_PathID=5
                unloadX=self.CuX
                unloadF=self.CuF
                Res=self.__General_unload(new_X,unloadX,unloadF)
            
        # on lower r4 path 3
        if PathID==3:
            tar_X=self.beta*self.tracker[3]
            if tar_X>-3*self.F1/self.k0:
                tar_X=-3*self.F1/self.k0
            
            if DX<0 and new_X<=tar_X:
                minX=new_X
                New_PathID=1
                Res=self.__Backbone(new_X)
            if DX<0 and new_X>tar_X:
                Res=self.__Pinch_lower(new_X)
            if DX>0:
                New_PathID=5
                unloadX=self.CuX
                unloadF=self.CuF
                Res=self.__General_unload(new_X,unloadX,unloadF)

        # on general unloading path 5
        if PathID==5:
            # calculate forces from path 2, 3, 5
            F_p2=self.__Pinch_upper(new_X)
            F_p3=self.__Pinch_lower(new_X)
            F_p5=self.CuF+(new_X-self.CuX)*self.k0*self.r3
            if F_p5<=F_p2 and F_p5>=F_p3:
                Res=F_p5
            if F_p5>F_p2:
                New_PathID=2
                Res=F_p2
            if F_p5<F_p3:
                New_PathID=3
                Res=F_p3
        
        # make sure there is no overshooting
        temp=self.__Backbone(new_X)
        if new_X>self.xu and Res>temp:
            Res=temp
        if new_X<-self.xu and Res<temp:
            Res=temp
        
        # print('NPath=',New_PathID,'NX=',new_X,'NY=',Res)
        # print('newK=',(Res-self.CuF)/(new_X-self.CuX))

        return Res
 
    
    def Estimate_tracker(self, new_X):
        maxX=self.tracker[0]
        unloadX_P=self.tracker[1]   #unloading form P backbone
        minX=self.tracker[2]
        unloadX_N=self.tracker[3]   #unloading from N backbone
        unloadX=self.tracker[4]     #normal unloading X
        unloadF=self.tracker[5]     #normal unloading F
        PathID=self.tracker[6]      #current location
        DX=new_X-self.CuX

        Res=0
        New_PathID=PathID  #if nothing trigers, then stay on current path

        # on positive back bone  0
        if PathID==0 and DX>=0:
            maxX=new_X
            Res=self.__Backbone(new_X)
        
        if PathID==0 and DX<0:
            unloadX_P=self.CuX
            New_PathID=4
            Res=self.__Backbone_unload_P(new_X,unloadX_P)
        
        # on negative back bone   1
        if PathID==1 and DX<=0:
            minX=new_X
            Res=self.__Backbone(new_X)
        
        if PathID==1 and DX>0:
            unloadX_N=self.CuX
            New_PathID=4
            Res=self.__Backbone_unload_N(new_X,unloadX_N)
        
        # on unloading path from backbone 4
        if PathID==4:
            if self.CuX>0:
                # find intersection of unloading from positive to pinching line
                inX,inY=Line_interXY(self.r3*self.k0,self.CuX,self.CuF,self.r4*self.k0,0,-self.F1)

                if new_X>unloadX_P:
                    maxX=new_X
                    New_PathID=0
                    Res=self.__Backbone(new_X)
                if new_X<unloadX_P and new_X>inX:
                    Res=self.__Backbone_unload_P(new_X,unloadX_P)
                if new_X<inX:
                    New_PathID=3
                    Res=self.__Pinch_lower(new_X)
            
            if self.CuX<0:
                # find intesection of unloading from negative to pinching line
                inX,inY=Line_interXY(self.r3*self.k0,self.CuX,self.CuF,self.r4*self.k0,0,self.F1)
                
                if new_X<unloadX_N:
                    minX=new_X
                    New_PathID=1
                    Res=self.__Backbone(new_X)
                if new_X>unloadX_N and new_X<inX:
                    Res=self.__Backbone_unload_N(new_X,unloadX_N)
                if new_X>inX:
                    New_PathID=2
                    Res=self.__Pinch_upper(new_X)

        # on upper r4 path 2
        if PathID==2:
            tar_X=self.beta*self.tracker[1]
            
            if DX>0 and new_X>tar_X:
                maxX=new_X
                New_PathID=0
                Res=self.__Backbone(new_X)
            if DX>0 and new_X<tar_X:
                Res=self.__Pinch_upper(new_X)
            if DX<0:
                New_PathID=5
                unloadX=self.CuX
                unloadF=self.CuF
                Res=self.__General_unload(new_X,unloadX,unloadF)
            
        # on lower r4 path 3
        if PathID==3:
            tar_X=self.beta*self.tracker[3]
            
            if DX<0 and new_X<tar_X:
                minX=new_X
                New_PathID=1
                Res=self.__Backbone(new_X)
            if DX<0 and new_X>tar_X:
                Res=self.__Pinch_lower(new_X)
            if DX>0:
                New_PathID=5
                unloadX=self.CuX
                unloadF=self.CuF
                Res=self.__General_unload(new_X,unloadX,unloadF)

        # on general unloading path 5
        if PathID==5:
            # calculate forces from path 2, 3, 5
            F_p2=self.__Pinch_upper(new_X)
            F_p3=self.__Pinch_lower(new_X)
            F_p5=self.CuF+(new_X-self.CuX)*self.k0*self.r3
            if F_p5<=F_p2 and F_p5>=F_p3:
                Res=F_p5
            if F_p5>F_p2:
                New_PathID=2
                Res=F_p2
            if F_p5<F_p3:
                New_PathID=3
                Res=F_p3
        
        # make sure there is no overshooting  # no need to bound force if just want to do tracker
        #temp=self.__Backbone(new_X)
        #if new_X>self.xu and Res>temp:
        #    Res=temp
        #if new_X<-self.xu and Res<temp:
        #    Res=temp

        return [maxX,unloadX_P,minX,unloadX_N,unloadX,unloadF,New_PathID]
    
    def GetK0(self):
        return self.k0
# Spring sub class Spr_EPHM   ID=4
class Spr_EPHM(Spring):
    def SetParameter(self, inputP):
        self.type=4
        self.parameter=inputP
        self.CuK=1  #set initial k as k0
        
    def init_tracker(self):  # there is no need to track
        self.tracker=0
        
    def GetNewForce(self,new_X):
        return 0
    
    def Estimate_tracker(self,new_X):
        return 0
    
    def GetK0(self):
        return self.parameter[0]
# Spring sub class Spr_MultiLinear    ID=5
class Spr_Multilinear(Spring):
    # Multilinear spring has N parameters
    # N_segment D1 F1 D2 F2 D3 F3....DN FN  (K beyond DN will be extrapolated)
    def SetParameter(self, inputP):
        self.type=5
        self.parameter=inputP
        self.CuK=self.parameter[1]/self.parameter[0]  #set initial k as k0
        
    def init_tracker(self):  # there is no need to track
        self.tracker=0
        
    def GetNewForce(self,new_X):
        temp=self.parameter[1:]
        X=temp[1::2]  # odd element
        F=temp[::2]     # even element
        temp=np.interp(abs(new_X),X,F)
        return temp*np.sign(new_X)
    
    def Estimate_tracker(self,new_X):
        return 0
    
    def GetK0(self):
        return self.parameter[0]
# Spring sub class Spr_CompOnly    ID=6
class Spr_CompOnly(Spring):
    # Comp only spring has 1 parameters
    # k0  (K beyond DN will be extrapolated)
    def SetParameter(self, inputP):
        self.type=6
        self.parameter=inputP
        self.CuK=self.parameter[0]/100000  #set initial k as k0/100000
        
    def init_tracker(self):  # there is no need to track
        self.tracker=0
        
    def GetNewForce(self,new_X):
        if new_X>0:
            return new_X*self.parameter[0]/100000
        else:
            return new_X*self.parameter[0]
    
    def Estimate_tracker(self,new_X):
        return 0
    
    def GetK0(self):
        return self.parameter[0]/100000
# Spring sub class Spr_TensionOnly   ID=7
class Spr_TensionOnly(Spring):
    # Tension only spring has 1 parameters
    # k0  (K beyond DN will be extrapolated)
    def SetParameter(self, inputP):
        self.type=7
        self.parameter=inputP
        self.CuK=self.parameter[0]  #set initial k as k0
        
    def init_tracker(self):  # there is no need to track
        self.tracker=0
        
    def GetNewForce(self,new_X):
        if new_X>0:
            return new_X*self.parameter[0]
        else:
            return new_X*self.parameter[0]/100000
    
    def Estimate_tracker(self,new_X):
        return 0
    
    def GetK0(self):
        return self.parameter[0]

# Model_file super class
class Model_file:
    def __init__(self):
        self.type=-1        # by default model type=-1
        pass        # create internal variables to store model info

    def LoadFile(self,lines):
        pass        # load from file

    def SaveFile(self,fileLoc):
        pass        # save to a file

    def To_str(self)-> str:
        pass

class Model_file_SDOF(Model_file):
    
    def __init__(self) -> None:
        self.type=1  # file type is 1 now
        self.mass=0
        self.spr_type=0
        self.spr_parameter=[]

    def LoadFile(self, lines):
        
        
        self.type=int(lines[0].strip())
        self.mass=float(lines[1].strip())
        self.spr_type=float(lines[2].strip())
        temp=lines[3].strip()
        self.spr_parameter=[float(num) for num in temp.split()]

        #SDOF file format
        #1 1 (indicating model type)
        #2 mass
        #3 Sprtype 1-linear  2-bilinear  3-CUREE
        #4 SprParameters 
        # 1 k0. 
        # 2 k0 ky Dy
        # 3 k0 xu F0 F1 r1 r2 r3 r4 alpha beta 

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
        tempstr+=str(self.type)+"\n"
        tempstr+=str(self.mass)+"\n"
        tempstr+=str(self.spr_type)+"\n"
        tempstr+="\t".join(str(f) for f in self.spr_parameter)
        return tempstr

class Model_file_sbld(Model_file):
    def __init__(self) -> None:
        self.type=2  # file type is 2 now
        self.mass=0
        self.spr_type=0
        self.spr_parameter=[]

    def LoadFile(self, lines):        
        
        self.type=int(lines[0].strip())
        temp=lines[1].strip()
        self.mass=[float(num) for num in temp.split()]
        temp=lines[2].strip()
        self.spr_type=[float(num) for num in temp.split()]
        
        for ii in range(len(self.mass)):
            temp=lines[ii+2].strip()
            self.spr_parameter[ii]=[float(num) for num in temp.split()]

        #sbld file format
        #1 2 (indicating model type)
        #2 mass1 mass2 mass3 ... massN
        #3 Sprtype1 Sprtype2 ... SprtypeN 1-linear  2-bilinear  3-CUREE
        #4 SprParameters[0] for sp1
        #5 SprParameters[1] for sp1
        #...
        #N+3 SprParameters[N-1] for spN

        # 1 k0. 
        # 2 k0 ky Dy
        # 3 k0 xu F0 F1 r1 r2 r3 r4 alpha beta 

    def SaveFile(self, fileLoc):
        with open(fileLoc, 'w') as file:
            file.write(self.To_str())
            
    def To_str(self) -> str:
        tempstr=""
        tempstr+=str(self.type)+"\n"
        for ii in range(len(self.mass)):
            tempstr+=str(self.mass[ii])
        tempstr+="\n"
        for ii in range(len(self.spr_type)):
            tempstr+=str(self.spr_type[ii])
        tempstr+="\n"

        for ii in range(len(self.spr_type)):
            tempstr+="\t".join(str(f) for f in self.spr_parameter[ii])
            tempstr+="\n"
        
        return tempstr
    
# Model_Dyn super class Model
class Model_Dyn:

    def __init__(self):
        self.type=-1
        pass        # create variable spaces

    def Construct(self,Modelfile:Model_file):
        pass        # construct variables based on model

    def Initialize(self):
        pass        # renew and make the model "new"

    def Analysis_NB(self, EQ:Earthquake,ScaleFactor,DampR, timeStep, Pro_Bar, WriteSpringOK):
        pass        # earthquake analysis

    def Analysis_Push(self,Pro:Protocols,DOF_ID,ScaleFactor):
        pass        # let's try this, it is a SDOF displacement-control push

    def DofPlot(self,DOF_ID,canvas):
        pass        # plot any DOF over time

    def HystPlot(self,Spr_ID,canvas):
        pass        # plot any spring element's hystersis

# Model sub class Model_SDOF
class Model_Dyn_SDOF(Model_Dyn):
    def __init__(self):
        self.type=1
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
            
            # turn model level X into spring level x
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

    def DofPlot(self, DOF_ID, canvas):
        
        X=self.time
        Y=self.GlobalX

        plot_on_canvas(X,Y,canvas)

        # print(len(X))
        # print(len(Y))
        # # Scale the data points to fit within the canvas
        # x_min, x_max = min(X), max(X)
        # y_min, y_max = min(Y), max(Y)
        # for i in range(len(X) - 1):
        #     # Map X and Y to canvas coordinates for each pair of adjacent points
        #     x1 = (X[i] - x_min) * canvas.winfo_width() / (x_max - x_min)
        #     y1 = canvas.winfo_height() - (Y[i] - y_min) * canvas.winfo_height() / (y_max - y_min)
        #     x2 = (X[i + 1] - x_min) * canvas.winfo_width() / (x_max - x_min)
        #     y2 = canvas.winfo_height() - (Y[i + 1] - y_min) * canvas.winfo_height() / (y_max - y_min)
        #     # Draw a line segment connecting adjacent data points
        #     canvas.create_line(x1, y1, x2, y2, fill="blue")

    def HystPlot(self, Spr_ID, canvas):
        
        self.Spr.HysPlot(canvas)  #SDOF model has only 1 spring, so Spr_ID is not used
   
# Model sub class Model_ShearBld
class Model_Dyn_sbld(Model_Dyn):
    def __init__(self):
        self.type=2
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

    def Construct(self, Modelfile: Model_file_sbld):
        self.mass=Modelfile.mass
        self.Spr_type=Modelfile.spr_type

        for ii in range(len(self.Spr_type)):
            self.Spr[ii]=Assign_Spr(self.Spr_type[ii])
            self.Spr[ii].SetParameter(Modelfile.spr_parameter[ii])

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
            
            # turn model level X into spring level x
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

    def DofPlot(self, DOF_ID, canvas):
        
        X=self.time
        Y=self.GlobalX

        plot_on_canvas(X,Y,canvas)

        # print(len(X))
        # print(len(Y))
        # # Scale the data points to fit within the canvas
        # x_min, x_max = min(X), max(X)
        # y_min, y_max = min(Y), max(Y)
        # for i in range(len(X) - 1):
        #     # Map X and Y to canvas coordinates for each pair of adjacent points
        #     x1 = (X[i] - x_min) * canvas.winfo_width() / (x_max - x_min)
        #     y1 = canvas.winfo_height() - (Y[i] - y_min) * canvas.winfo_height() / (y_max - y_min)
        #     x2 = (X[i + 1] - x_min) * canvas.winfo_width() / (x_max - x_min)
        #     y2 = canvas.winfo_height() - (Y[i + 1] - y_min) * canvas.winfo_height() / (y_max - y_min)
        #     # Draw a line segment connecting adjacent data points
        #     canvas.create_line(x1, y1, x2, y2, fill="blue")

    def HystPlot(self, Spr_ID, canvas):
        
        self.Spr.HysPlot(canvas)  #SDOF model has only 1 spring, so Spr_ID is not used
   
# Model sub class Model_Biaxial

# Model sub class Model_6DOF
