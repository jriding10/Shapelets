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
import time
import matplotlib
import numpy as np
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

options = dict([('title',None), ('dirname', None), ('filetypes', None), ('fits_data', None),
                ('nmax', 5), ('im_array', None), ('sizeData', None)])
drawing = dict([('x0', 0), ('y0', 0), ('x1', None), ('y1', None)])

def quitter():
    global root
    if messagebox.askyesno('Verify', 'Do you really want to quit?'):
        plt.close('all')
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
    basis_entry.config(fg='black')
    root.update_idletasks()

def getFitsFile():
    fits_data = pyfits.open(options['fits_data'])
    source = fits_data[0].data
    time.sleep(1)
    fits_data.close()
    checkSourceData(source)
    displayFitsData()

def checkSourceData(source):
    global file_label, root
    size_of_data = source.shape
    dim_of_data = len(size_of_data)

    if (dim_of_data == 2):
        options['im_array'] = source
        options['sizeData'] = options['im_array'].shape
    elif (dim_of_data < 2):
        file_label.config(text='Invalid fits file')
        file_label.config(bg = 'red')
        root.update_idletasks()
        return
    else:
        options['im_array'] = source[0,0,:,:]
        options['sizeData'] = options['im_array'].shape

def xy(event):
    global lastx, lasty
    lastx, lasty = event.xdata, event.ydata
    drawing['x0'] = lastx
    drawing['y0'] = lasty

def selectROI(event):
    global lastx, lasty
    x, y = event.xdata, event.ydata
    drawing['x1'] = x
    drawing['y1'] = y
    deltax = x - lastx
    deltay = y - lasty
    rect = plt.Rectangle((lastx,lasty),deltax,deltay,linewidth=2,edgecolor='k', facecolor='none')
    ax.add_patch(rect)
    fig.canvas.draw()

def drawRectangle():
    x = 0
    y = 0
    rect = plt.Rectangle((x,y), width=100, height=100, edgecolor='g')
    ax.add_patch(rect)
    fig.canvas.draw()


def displayFitsData():
    ax.clear()
    dataset = options['im_array']
    plt.imshow(dataset, cmap='hot')
    fig.canvas.draw()

def updateFitsDisplay():
    ax.clear()
    xextent = int(drawing['x1'] - drawing['x0'])
    yextent = int(drawing['y1'] - drawing['y0'])
    new_source = np.zeros((yextent, xextent))

    for i in range(yextent):
        for j in range(xextent):
            rows = int(i+int(drawing['y0']))
            cols = int(j+int(drawing['x0']))
            new_source[i, j] = options['im_array'][rows, cols]

    plt.imshow(new_source, 'hot')
    fig.canvas.draw()

########################################################################################################################
root = Tk()
root.title('Playing with GUIs')

lastx, lasty = 0, 0

set_fonts = ('Helvetica', 12, 'bold')

file_label = Label(root, fg='black', justify='left', font='Helvetica 12')
get_file = Button(root, text='Select', command=getFilename)
basis_label = Label(root, text='Number of basis:', fg='black', font='Helvetica 12')
basis_entry = Entry(root, width=3, fg='red')
basis_entry.bind('<Return>', getNumberBasis)
quitter = Button(root, text='Quit', command=quitter)
load_file = Button(root, text='Load', command=getFitsFile)
draw_rect = Button(root, text='Draw', command=drawRectangle)
resize_im = Button(root, text='Resize', command=updateFitsDisplay)
undo_resize = Button(root, text='Undo', command=displayFitsData)

fig = plt.figure(1)
ax = fig.add_subplot(111)
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(column=1, row=3, sticky=(N,W,E,S))

idpress = fig.canvas.mpl_connect('button_press_event', xy)
idmove = fig.canvas.mpl_connect('button_release_event', selectROI)
# #idrel = fig.canvas.mpl_connect('button_release_event', )

get_file.grid(column=0, row=0)
file_label.grid(column=1, row=0, sticky='W')
load_file.grid(column=2, row=0)
basis_label.grid(column=0, row=1, sticky='E')
basis_entry.grid(column=1, row=1, sticky='W')
quitter.grid(column=0, row=2)
resize_im.grid(column=2, row=1)
draw_rect.grid(column=1, row=2)
undo_resize.grid(column=2, row=2)


mainloop()