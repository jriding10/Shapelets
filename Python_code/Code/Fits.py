import astropy.io.fits as pyfits
import numpy as np
import sys

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
        self.yAngScale = np.radians(fitsHeader['CDELT1'])
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
            sys.exit()
