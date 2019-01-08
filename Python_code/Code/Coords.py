import numpy as np

class Coordinates:
    def __init__(self, colRA_pxl=None,  rowDec_pxl=None, xextent=None, yextent=None, pxl_coords=None,
                 raAxis=[], decAxis=[]):
        self.colRA_pxl = colRA_pxl
        self.rowDec_pxl = rowDec_pxl
        self.xextent = xextent
        self.yextent = yextent
        self.pxl_coords = pxl_coords
        self.raAxis = raAxis
        self.decAxis = decAxis

    def changeRADecRef(self, x, position, pxl_position, angScale):
        new_row = x[0]
        new_column = x[1]
        deltaRowPixel = abs(pxl_position[0] - new_row)
        deltaColumnPixel = abs(pxl_position[1] - new_column)
        newRA = position[0] - angScale[1]*deltaColumnPixel
        newDec = position[1] - angScale[0]*deltaRowPixel
        self.colRA_pxl = new_column
        self.rowDec_pxl = new_row
        position[0] = newRA
        position[1] = newDec
        return position

    def calcCoordinateList(self, angScale):
        total_coordinates = self.xextent * self.yextent
        coordinates = np.zeros((total_coordinates, 2))
        row_start = int(-1 * self.xextent/2)
        column_start = int(-1 * self.yextent/2)

        k = 0
        for i in range(self.xextent):
            for j in range(self.yextent):
                coordinates[k, 0] = (row_start + i) * abs(angScale[1])
                coordinates[k, 1] = (column_start + j) * abs(angScale[0])
                k += 1
        self.pxl_coords = coordinates
        return coordinates

    def expandArray(self, flat_array):
        num_rows = self.xextent
        num_cols = self.yextent
        new_array = np.zeros((num_rows, num_cols))
        k = 0

        for i in range(num_rows):
            for j in range(num_cols):
                new_array[i,j] = flat_array[k]
                k += 1

        return new_array

    def resizeData(self, source, x1, x2):
       self.xextent = int(x2[0] - x1[0])
       self.yextent = int(x2[1] - x1[1])
       temp_source = np.zeros((self.xextent, self.yextent))

       for i in range(self.xextent):
           for j in range(self.yextent):
               temp_source[i,j] = source[i+x1[0], j+x1[1]]

       return temp_source

    def raCoords(self, position, angScale):
        startRA = position[0] + self.colRA_pxl*angScale[0]
        self.raAxis = np.zeros((self.yextent))
        for i in range(self.yextent):
            self.raAxis[i] = np.degrees(startRA - i*angScale[0])
        return self.raAxis

    def decCoords(self, position, angScale):
        startDec = position[1] - self.rowDec_pxl*angScale[1]
        self.decAxis = np.zeros((self.xextent))
        for i in range(self.xextent):
            self.decAxis[i] = np.degrees(startDec + i*angScale[1])
        return self.decAxis
