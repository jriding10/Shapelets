import numpy as np

class Coordinates:
    def __init__(self, colRA_pxl=None,  rowDec_pxl=None, xextent=None, yextent=None,
                 raAxis=[], decAxis=[], data = []):
        self.colRA_pxl = colRA_pxl
        self.rowDec_pxl = rowDec_pxl
        self.xextent = xextent
        self.yextent = yextent
        self.raAxis = raAxis
        self.decAxis = decAxis
        self.data = data

    def changeRADecRef(self, old_position, new_position, ra_decPosition, angScale):
        new_row = new_position[0]
        new_column = new_position[1]
        deltaRowPixel = abs(old_position[0] - new_row)
        deltaColumnPixel = abs(old_position[1] - new_column)
        newRA = ra_decPosition[0] - angScale[1]*deltaColumnPixel
        newDec = ra_decPosition[1] - angScale[0]*deltaRowPixel
        self.colRA_pxl = new_column
        self.rowDec_pxl = new_row
        ra_decPosition[0] = newRA
        ra_decPosition[1] = newDec
        return ra_decPosition

    def calcTransformCoordinates(self, pa, coords):
        transform_matrix = [[m.cos(pa), -m.sin(pa)], [m.sin(pa), m.cos(pa)]]
        coord_list = np.dot(transform_matrix, np.transpose(coords))
        coords = np.transpose(coord_list)
        return coords

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

    def useWholeImage(self, extent, dataset, old_pxls, old_coords, angScale):
        self.row0 = 0
        self.column0 = 0
        self.row1 = extent[1]
        self.column1 = extent[0]
        self.row_extent = abs(self.row1 - self.row0)
        self.column_extent = abs(self.column1 - self.column0)
        ref_coords = self.selectROI(dataset, old_pxls, old_coords, angScale)
        return ref_coords

    def selectROI(self, dataset, old_pxls, old_coords, angScale):
        new_source = np.zeros((self.row_extent, self.column_extent))
        new_ref_pxl = [row0, col0]
        ref_coords = self.changeRADecRef(old_pxls, new_ref_pxl, old_coords, angScale)

        for i in range(self.row_extent):
            for j in range(self.column_extent):
                rows = int(i + self.row0)
                cols = int(j + self.column0)
                new_source[i, j] = dataset[rows, cols]

        self.colRA_pxl = 0
        self.rowDec_pxl = 0
        self.data = new_source
        return ref_coords

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
