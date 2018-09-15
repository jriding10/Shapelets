#! /bin/bash python

# so I can play with gui functions before using them with the real code
# need: file dialog and then display that filepath
#       display/canvas
#       way to specify number of basis
#       way to show goodness of fit
#       way to save text and jpgs
#       quit
#       resize and undo
#       go

import os
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
import astropy.io.fits as pyfits
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure

options = dict([('title',None), ('dirname', None), ('filetypes', None), ('fits_data', None),
                ('nmax', 5), ('im_array', None)])

def quitter():
    global root
    if messagebox.askyesno('Verify', 'Do you really want to quit?'):
        root.destroy()

def getFilename():
    options['dirname'] = os.getcwd()
    options['title'] = 'Select fits file'
    options['filetypes'] = (('fits files', '*.fits'), ('all files', '*.*'))

    fits_file = filedialog.askopenfilename(initialdir=options['dirname'], title=options['title'],
                                          filetypes=options['filetypes'])
    options['fits_data'] = fits_file
    updateDisplay()
    return

def updateDisplay():
    global root
    global file_label
    file_label.config(text=options['fits_data'])
    root.update_idletasks()

def getNumberBasis(event):
    global basis_entry
    global root
    n = basis_entry.get()
    options['nmax'] = int(n)
    basis_entry.config(fg='green')
    root.update_idletasks()

def getFitsFile():
    fits_data = pyfits.open(options['fits_data'])
    source = fits_data[0].data
    options['im_array'] = source
    fits_data.close()
    return

def displayFitsData():
    dataset = options['im_array']
    plt.contour(dataset)


########################################################################################################################
root = Tk()
root.title('Playing with GUIs')

set_fonts = ('Helvetica', 12, 'bold')

file_label = Label(root, fg='black', justify='left', font='Helvetica 12', width=50)
get_file = Button(root, text='Select', command=getFilename)
basis_label = Label(root, text='Number of basis:', fg='black', font='Helvetica 12')
basis_entry = Entry(root, width=3)
basis_entry.bind('<Return>', getNumberBasis)
quitter = Button(root, text='Quit', command=quitter)
load_file = Button(root, text='Load', command=getFitsFile)

fits_image = Figure(figsize=(5, 4), dpi=100)
graph_handle = fits_image.add_subplot(111)
dataset = options['im_array']
plt.contour(dataset)

canvas = FigureCanvasTkAgg(fits_image, master=root)
canvas.show()
canvas.get_tk_widget().grid(column=0, row=3)

toolbar = NavigationToolbar2TkAgg(canvas, root)
toolbar.update()
canvas._tkcanvas.grid(column=0, row=4)


get_file.grid(column=0, row=0)
file_label.grid(column=1, row=0, sticky='W')
load_file.grid(column=2, row=0)
basis_label.grid(column=0, row=1, sticky='E')
basis_entry.grid(column=1, row=1, sticky='W')
quitter.grid(column=0, row=2)


mainloop()