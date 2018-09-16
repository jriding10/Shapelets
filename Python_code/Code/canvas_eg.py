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
from tkinter import ttk
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
    basis_entry.config(fg='black')
    root.update_idletasks()

def getFitsFile():
    fits_data = pyfits.open(options['fits_data'])
    source = fits_data[0].data
    options['im_array'] = source
    fits_data.close()

def xy(event):
    global lastx, lasty
    lastx, lasty = event.x, event.y

def addLine(event):
    global lastx, lasty
    x, y = canvas.canvasx(event.x), canvas.canvasy(event.y)
    canvas.create_line((lastx, lasty, x, y), width=5, fill=color, tags='currentline')
    lastx, lasty = x, y

# see eg location 1622 of Modern tkinter
def setColor(newcolor):
    global color
    color = newcolor
    canvas.dtag('all', 'paletteSelected')
    canvas.itemconfigure('palette', outline='white')
    canvas.addtag('paletteSelected', 'withtag', 'palette%s' % color)
    canvas.itemconfigure('paletteSelected', outline='#999999')

def doneStroke(event):
    canvas.itemconfigure('currentline', width=1)

########################################################################################################################
root = Tk()
root.title('Playing with GUIs')

lastx, lasty = 0, 0
color = 'black'


set_fonts = ('Helvetica', 12, 'bold')

file_label = Label(root, fg='black', justify='left', font='Helvetica 12', width=50)
get_file = Button(root, text='Select', command=getFilename)
basis_label = Label(root, text='Number of basis:', fg='black', font='Helvetica 12')
basis_entry = Entry(root, width=3, fg='red')
basis_entry.bind('<Return>', getNumberBasis)
quitter = Button(root, text='Quit', command=quitter)
load_file = Button(root, text='Load', command=getFitsFile)

h=ttk.Scrollbar(root, orient=HORIZONTAL)
v=ttk.Scrollbar(root, orient=VERTICAL)

canvas = Canvas(root, scrollregion=(0,0,1000,1000), yscrollcommand=v.set, xscrollcommand=h.set)
h['command'] = canvas.xview
v['command'] = canvas.yview
ttk.Sizegrip(root).grid(column=2, row=4, sticky=(S,E))

canvas.grid(column=1, row=3, sticky=(N,W,E,S))
h.grid(column=1, row=4, sticky=(W,E))
v.grid(column=2, row=3, sticky=(N,S))
root.grid_columnconfigure(0, weight=1)
root.grid_rowconfigure(0, weight=1)

canvas.bind("<Button-1>", xy)
canvas.bind("<B1-Motion>", addLine)
canvas.bind('<B1-ButtonRelease>', doneStroke)

id = canvas.create_rectangle((10,10,30,30), fill='red', tags=('palette', 'plattered'))
canvas.tag_bind(id, '<Button-1>', lambda x: setColor('red'))
id = canvas.create_rectangle((10,35,30,55), fill='blue', tags=('palette', 'platteblue'))
canvas.tag_bind(id, '<Button-1>', lambda x: setColor('blue'))
id = canvas.create_rectangle((10,60,30,80), fill='black', tags=('palette', 'platteblack', 'paletteSelected'))
canvas.tag_bind(id, '<Button-1>', lambda x: setColor('black'))

setColor('black')
canvas.itemconfigure('palette', width=5)


get_file.grid(column=0, row=0)
file_label.grid(column=1, row=0, sticky='W')
load_file.grid(column=2, row=0)
basis_label.grid(column=0, row=1, sticky='E')
basis_entry.grid(column=1, row=1, sticky='W')
quitter.grid(column=0, row=2)


mainloop()