#! /bin/bash/python

# creates pretty pictures
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np


class FitsData:
    def __init__(self, filename=None, RA=None, Dec=None, xAngScale=None, yAngScale=None, rowDecpxl=None,
                 colRApxl=None, xExtent=None, yExtent=None, source_data=None, raAxis=None, decAxis=None):
        self.filename = filename
        self.RA = RA
        self.Dec = Dec
        self.xAngScale = xAngScale
        self.yAngScale = yAngScale
        self.rowDecpxl = rowDecpxl
        self.colRApxl = colRApxl
        self.xExtent = xExtent
        self.yExtent = yExtent
        self.source_data = source_data
        self.raAxis = raAxis
        self.decAxis = decAxis

    def getFitsFile(self):
        fitsFile = pyfits.open(self.filename)
        fitsHeader = fitsFile[0].header
        fitsData = fitsFile[0].data

        # RA == y axis/column indices, Dec == x axis/row indices
        self.RA = np.radians(fitsHeader['CRVAL1'])
        self.Dec = np.radians(fitsHeader['CRVAL2'])
        self.xAngScale = np.radians(fitsHeader['CDELT2'])
        self.yAngScale = np.radians(fitsHeader['CDELT2'])
        self.rowDecpxl = fitsHeader['CRPIX2']
        self.colRApxl = fitsHeader['CRPIX1']
        self.xExtent = fitsHeader['NAXIS2']
        self.yExtent = fitsHeader['NAXIS1']
        fitsFile.close()

        self.checkSourceData(fitsData)

    def checkSourceData(self, fitsdata):
        size_of_data = fitsdata.shape
        dim_of_data = len(size_of_data)

        if dim_of_data == 2:
            self.source_data = fitsdata
        elif dim_of_data == 4:
            self.source_data = fitsdata[0, 0, :, :]
        else:
            print('Fits data is of unknown dimensions: %.0f' % dim_of_data)

    def resizeData(self):
        self.changeRADecRef(x1, y1)
        self.xExtent = int(x2 - x1)
        self.yExtent = int(y2 - y1)
        temp_source = np.zeros((self.xExtent, self.yExtent))

        for i in range(self.xExtent):
            for j in range(self.yExtent):
                temp_source[i,j] = self.source_data[i+x1, j+y1]

        self.source_data = temp_source

    def raCoords(self):
        startRA = self.RA + self.colRApxl*self.yAngScale
        self.raAxis = np.zeros((self.yExtent))
        for i in range(self.yExtent):
            self.raAxis[i] = np.degrees(startRA - i*self.yAngScale)

    def decCoords(self):
        startDec = self.Dec - self.rowDecpxl*self.xAngScale
        self.decAxis = np.zeros((self.xExtent))
        for i in range(self.xExtent):
            self.decAxis[i] = np.degrees(startDec + i*self.xAngScale)

    def changeRADecRef(self, new_row, new_column):
        new_row = int(new_row)
        new_column = int(new_column)
        deltaRowPixel = abs(self.rowDecpxl - new_row)
        deltaColumnPixel = abs(self.colRApxl - new_column)
        self.RA = self.RA - self.yAngScale*deltaColumnPixel
        self.Dec = self.Dec - self.xAngScale*deltaRowPixel
        self.colRApxl = 0
        self.rowDecpxl = 0

fits = FitsData()
fits.filename = '../Fits_files/FnxA.fits'
fits.getFitsFile()

#3C353
x1, y1 = 100, 100
x2, y2 = 250, 250
ramin = 259.798
ramax = 260.005
decmin = -1.092
decmax = -0.885

#PKS0410-75
x1, y1 = 150, 150
x2, y2 = 210, 210
ramin = 62.077
ramax = 62.159
decmin = -75.165
decmax = -75.083

#Fornax A
x1, y1 = 275, 150
x2, y2 = 825, 900
ramin = 49.126
ramax = 50.167
decmin = -37.539
decmax = -36.776

fits.resizeData()
fits.raCoords()
fits.decCoords()

plt.figure(0)
plt.imshow(fits.source_data, cmap='hot', extent=[ramax, ramin, decmin, decmax])
plt.xlabel('RA (degrees)')
plt.ylabel('Dec (degrees)')
plt.colorbar()
plt.show()



