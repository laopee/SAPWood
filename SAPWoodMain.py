# This is the main program including GUI
# this contains sub-routines to read EQ and load files into obj
# sub-routines to load model files into model obj
# sub routines to do push over or time history analysis or time history analysis
# GUIs to control analysis and display results
import SAPWoodClass as SP

import tkinter as tk
import matplotlib.pyplot as plt

#Global variables
EQ_current=SP.Earthquake()



#Button functions
def load_eq():
    filename=tk.filedialog.askopenfilename()
    EQ_current.LoadEQ(filename)
    print(EQ_current.Ax[3])
    print(EQ_current.Ay[5])

def PlotAx():
    plt.plot(EQ_current.t,EQ_current.Ax)
    plt.show()

def PlotAy():
    plt.plot(EQ_current.t,EQ_current.Ay)
    plt.show()

#All functions defined above, now do the GUIs

#GUIs
root=tk.Tk()
root.title=("load EQ file")

button_load=tk.Button(root,text="load", command=load_eq)
button_load.pack()

button_Px=tk.Button(root,text='plotX',command=PlotAx)
button_Py=tk.Button(root,text='plotY',command=PlotAy)
button_Px.pack()
button_Py.pack()

root.mainloop()