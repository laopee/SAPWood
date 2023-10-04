# This is the main program including GUI
# this contains sub-routines to read EQ and load files into obj
# sub-routines to load model files into model obj
# sub routines to do push over or time history analysis or time history analysis
# GUIs to control analysis and display results
import SAPWoodClass as SP

import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

#Global variables
EQ_current=SP.Earthquake()
Pro_current=SP.Protocols()



#Button functions
def load_eq():
    filename=tk.filedialog.askopenfilename()
    EQ_current.LoadEQ(filename)
    print(EQ_current.Ax[3])
    print(EQ_current.Ay[5])

def PlotAx():
    plot_xy_on_existing_canvas(Figure_Ax,EQ_current.t,EQ_current.Ax)

def PlotAy():
    plt.plot(EQ_current.t,EQ_current.Ay)
    plt.show()

#utility functions
def plot_xy_on_existing_canvas(canvas, X, Y):
    """
    Plots the X-Y data points as a line plot on an existing Tkinter canvas.

    Args:
        canvas (tk.Canvas): The existing canvas where the plot will be drawn.
        X (list): List of X coordinates.
        Y (list): List of corresponding Y coordinates.
    """
    # Scale the data points to fit within the canvas
    x_min, x_max = min(X), max(X)
    y_min, y_max = min(Y), max(Y)

    for i in range(len(X) - 1):
        # Map X and Y to canvas coordinates for each pair of adjacent points
        x1 = (X[i] - x_min) * canvas.winfo_width() / (x_max - x_min)
        y1 = canvas.winfo_height() - (Y[i] - y_min) * canvas.winfo_height() / (y_max - y_min)
        x2 = (X[i + 1] - x_min) * canvas.winfo_width() / (x_max - x_min)
        y2 = canvas.winfo_height() - (Y[i + 1] - y_min) * canvas.winfo_height() / (y_max - y_min)

        # Draw a line segment connecting adjacent data points
        canvas.create_line(x1, y1, x2, y2, fill="blue")



#All functions defined above, now do the GUIs

#GUIs
root=tk.Tk()
root.title=("NB Analysis")

# create all Tabs you need
tabControl=ttk.Notebook(root)
tab1=ttk.Frame(tabControl)
tab2=ttk.Frame(tabControl)
tab3=ttk.Frame(tabControl)
tab4=ttk.Frame(tabControl)

tabControl.add(tab1,text='Analysis')
tabControl.add(tab2,text='Results')
tabControl.add(tab3,text='Funny stuff')
tabControl.add(tab4,text='Push spring')

tabControl.pack(expand=1, fill="both")

# add controls in each tab
# tab 1
button_load=tk.Button(tab1,text="load", command=load_eq)
button_load.pack()

button_Px=tk.Button(tab1,text='plotX',command=PlotAx)
button_Py=tk.Button(tab1,text='plotY',command=PlotAy)
button_Px.pack()
button_Py.pack()

Figure_Ax = tk.Canvas(tab1, width=400, height=300)
Figure_Ax.pack(fill=tk.X)

# tab 2
button_save=tk.Button(tab2,text="Save")
button_save.pack()
# tab 3
button_game=tk.Button(tab3,text="Random ADV")
button_game.pack()

# tab 4
button_push=tk.Button(tab4,text="push it")
button_push.pack()

text4 = tk.Text(tab4, font=("Purisa", 12),height=4,width=40)
text4.pack()

Figure_Hys = tk.Canvas(tab4, width=400, height=300)
Figure_Hys.pack(fill=tk.X)

root.mainloop()