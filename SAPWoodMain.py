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
    plt.plot(EQ_current.t,EQ_current.Ax)
    plt.show()

#GUIs
root=tk.Tk()
root.title=("load EQ file")

button_load=tk.Button(root,text="load", command=load_eq)
button_load.pack()

root.mainloop()