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
Model_current=SP.Model_Dyn()
ModelFile_current=SP.Model_file()
Spr_current=SP.Spring()

#Global functions



# Assign model based on type from model file  #I would like to get these in with SAPWood Class, but I have trouble import functions from that
def Assign_Model(type:int):
    if type==1:
        return SP.Model_Dyn_SDOF()
    if type==2:
        return SP.Model_Dyn()    
#________________above functions are not ideal, may improve later

#Button functions
def load_eq():
    filename=tk.filedialog.askopenfilename()
    EQ_current.LoadEQ(filename)
    #print(EQ_current.Ax[3])
    #print(EQ_current.Ay[5])

def PlotAx():
    plot_xy_on_existing_canvas(Figure_Ax,EQ_current.t,EQ_current.Ax)

def PlotAy():
    plt.plot(EQ_current.t,EQ_current.Ay)
    plt.show()

def load_Mfile():
    global ModelFile_current
    global Model_current

    filename=tk.filedialog.askopenfilename()
    with open(filename,'r') as file:
            lines=file.readlines()
    # determine what model file type it is
    mytype=int(lines[0].strip())
    if mytype==1:
        ModelFile_current=SP.Model_file_SDOF()
    if mytype==2:
        ModelFile_current=SP.Model_file()


    ModelFile_current.LoadFile(lines)
    print('file loaded'+str(ModelFile_current.type))
    print(ModelFile_current.To_str())
    Model_current=Assign_Model(ModelFile_current.type) # this supposed to assign model type by input file

    Model_current.Construct(ModelFile_current)

def show_Mfile():
    msg=ModelFile_current.To_str()
    
    display_msg_on_existing_text(Tb1_text,msg)

def Ansys_NB():
    global Model_current, EQ_current

    Model_current.Analysis_NB(EQ_current,300,0.05,0.01,progress_bar,True)

def Plot_X():
    global Model_current

    Model_current.DofPlot(1,Figure_Ax)  
    print('I plotted X')  

def Plot_hys():
    global Model_current

    Model_current.HystPlot(1,Figure_Hys)
    print('now Hyst plot')
    print(max(Model_current.GlobalX))

def Spring_Push():
    global Spr_current, Pro_current
    # 

#__________end of button functions

#utility functions
def plot_xy_on_existing_canvas(canvas:tk.Canvas, X, Y):
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

def display_msg_on_existing_text(mytext:tk.Text,msg:str):
    mytext.delete('1.0','end')
    mytext.insert('1.0',msg)
    

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

button_loadM=tk.Button(tab1,text="load Model",command=load_Mfile)
button_loadM.pack()

Tb1_text=tk.Text(tab1,height=8)
Tb1_text.pack()

button_showModel=tk.Button(tab1,text="show model",command=show_Mfile)
button_showModel.pack()

Figure_Ax = tk.Canvas(tab1, width=200, height=100)
Figure_Ax.pack(fill=tk.X)

# analysis gadgets
values = ['0.05', '0.02', '0.01']
drop_down = ttk.Combobox(tab1, values=values)
drop_down.pack()

progress_bar = ttk.Progressbar(tab1, orient='horizontal', length=200, mode='determinate')
progress_bar.pack()

button_anlysis=tk.Button(tab1,text="Anlysis", command=Ansys_NB)
button_anlysis.pack()

button_response=tk.Button(tab1,text="X-result",command=Plot_X)
button_response.pack()

button_hyst=tk.Button(tab1,text='plot hys', command=Plot_hys)
button_hyst.pack()

# tab 2
button_save=tk.Button(tab2,text="Save")
button_save.pack()
# tab 3
button_game=tk.Button(tab3,text="Random ADV")
button_game.pack()

# tab 4
button_push=tk.Button(tab4,text="push it",command=Plot_hys)
button_push.pack()

text4 = tk.Text(tab4, font=("Purisa", 12),height=4,width=40)
text4.pack()

Figure_Hys = tk.Canvas(tab4, width=400, height=300)
Figure_Hys.pack(fill=tk.X)

root.mainloop()