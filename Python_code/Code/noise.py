#! /bin/bash python
# Noise.py

##########################################################################
## Does the full shapelet build etc but without the gui interface. The  ##
## idea is to examine the results in terms of noise performance.        ##
##########################################################################

import os
import math as m
import numpy as np
import matplotlib.pyplot as plt

import sys
import Fits
import Source
import Coords
import Shapes


# Pick a source: PKS0410, 3C353, ForA
which_source = 'PKS0410'
plot_zoomed_image = False
plot_image = False
plot_noise = False
plot_noise_hist = False
do_circular_fit = False
do_model_noise_cf = False
do_elliptical_fit = False
do_model_noise_cf2 = False
plot_circular_moments = False
plot_elliptical_moments = False
plot_nmse_circular = True
plot_nmse_elliptical = False


class Common:
    def __init__(self):
        self.nmax = 10
        self.x1 = [0, 0]
        self.x2 = [0, 0]
        self.n1 = [0, 0]
        self.n2 = [0, 0]
        self.position = [0., 0.]
        self.pxl_position = [0, 0]
        self.angScale = [0., 0.]
        self.coord_list = []
        self.raAxis = []
        self.decAxis = []
        self.dataset = []
        self.debug = None
        self.debug2 = None
        self.debug3 = None


def pretty_plot(dataset):
    ra_max = max(data.raAxis)
    ra_min = min(data.raAxis)
    dec_max = max(data.decAxis)
    dec_min = min(data.decAxis)

    plt.figure(0)
    plt.imshow(dataset, cmap='hot', extent=[ra_max, ra_min, dec_min, dec_max])
    plt.xlabel('RA (degrees)')
    plt.ylabel('Dec (degrees)')
    plt.colorbar()
    plt.show()

def moment_map(moments):
    plt.imshow(moments, cmap='coolwarm')
    plt.xlabel('n1')
    plt.ylabel('n2')
    plt.colorbar()
    plt.show()

def zoomed_image():
    data.position = coords.changeRADecRef(data.x1, data.position, data.pxl_position, data.angScale)
    data.pxl_position = data.x1
    dataset = coords.resizeData(data.dataset, data.x1, data.x2)
    data.coord_list = coords.calcCoordinateList(data.angScale)
    data.raAxis = coords.raCoords(data.position, data.angScale)
    data.decAxis = coords.decCoords(data.position, data.angScale)
    return dataset

def straight_image():
    data.x1 = [0, 0]
    coords.xextent = fits.xExtent
    coords.yextent = fits.yExtent
    data.position = coords.changeRADecRef(data.x1, data.position, data.pxl_position, data.angScale)
    data.pxl_position = data.x1
    data.coord_list = coords.calcCoordinateList(data.angScale)
    data.raAxis = coords.raCoords(data.position, data.angScale)
    data.decAxis = coords.decCoords(data.position, data.angScale)

def noise_image(dataset):
    data.position = coords.changeRADecRef(data.n1, data.position, data.pxl_position, data.angScale)
    data.pxl_position = data.n1
    noise = coords.resizeData(dataset, data.n1, data.n2)
    data.coord_list = coords.calcCoordinateList(data.angScale)
    data.raAxis = coords.raCoords(data.position, data.angScale)
    data.decAxis = coords.decCoords(data.position, data.angScale)
    return noise

def histogram_plot(dataset, nbins):
    dataset = dataset.flatten()
    plt.hist(dataset, bins=nbins)
    plt.xlabel('Value')
    plt.ylabel('Count')
    plt.show()

def histogram_2plot(dataset, model, nbins):
    dataset = dataset.flatten()
    model = model.flatten()
    data_max = max(dataset)
    data_min = min(dataset)
    binning = np.linspace(data_min, data_max, nbins)
    plt.hist(dataset, bins=binning, alpha=0.5, label='Noise')
    plt.hist(model, bins=binning, alpha=0.5, label='Model')
    plt.xlabel('Value')
    plt.ylabel('Count')
    plt.legend()
    plt.show()

def circular_fit(dataset):
    extent = [0, 0]
    extent[0] = int(data.x2[0] - data.x1[0])
    extent[1] = int(data.x2[1] - data.x1[0])
    major, minor, pa = source.calcDefaultValues(data.angScale, extent)
    shapes.major = minor
    shapes.minor = minor
    shapes.position_angle = 0.0
    model, resids = create_shapelet_fit(dataset)
    return model, resids

def elliptical_fit(dataset):
    extent = [0, 0]
    extent[0] = int(data.x2[0] - data.x1[0])
    extent[1] = int(data.x2[1] - data.x1[0])
    major, minor, pa = source.calcDefaultValues(data.angScale, extent)
    shapes.major = major
    shapes.minor = minor
    shapes.position_angle = pa
    model, resids = create_shapelet_fit(dataset)
    return model, resids

def create_shapelet_fit(dataset):
    shapes.calcMoments(dataset, data.coord_list, data.nmax)
    flat_model = shapes.calcShapeletModel(data.coord_list, shapes.moments)
    flat_data = dataset.flatten()
    flat_resids = shapes.calcResidual(flat_data, flat_model)
    model = coords.expandArray(flat_model)
    resids = coords.expandArray(flat_resids)
    return model, resids

def add_moments(dataset, noise):
    flat_data = dataset.flatten()
    data.debug2 = flat_data
    flat_noise = noise.flatten()
    noise_var = np.var(flat_noise)
    total_moments = shapes.moments.shape[0]
    nmse = np.zeros((total_moments, 2))
    ssim = np.zeros((total_moments, 1))
    for i in range(total_moments):
        moments = shapes.minimiseMoments(i+1)
        data.debug3 = moments
        flat_model = shapes.calcShapeletModel(data.coord_list, moments)
        flat_resids = shapes.calcResidual(flat_data, flat_model)
        nmse[i, 0] = source.calcNMSE(flat_data, flat_resids)
        nmse[i, 1] = noise_var
        ssim[i] = source.calcSSIM(flat_data, flat_model)

    plt.plot(nmse[:, 0], 'r', label='NMSE')
    plt.plot(nmse[:, 1], 'b', label='Variance of Noise')
    plt.legend()
    plt.show()
    plt.plot(ssim)
    plt.show()


if __name__ == '__main__':
    data = Common()
    fits = Fits.FitsData()
    shapes = Shapes.Shapelets()
    coords = Coords.Coordinates()
    source = Source.TheSource()

    if which_source == 'PKS0410':
        #PKS0410-75
        filename = '../Fits_files/PKS0410-75.fits'
        data.x1 = [150, 150]
        data.x2 = [210, 210]
        data.n1 = [100, 100]
        data.n2 = [150, 150]
    elif which_source == '3C353':
        #3C353
        filename = '../Fits_files/3C353.fits'
        data.x1 = [100, 100]
        data.x2 = [250, 250]
        data.n1 = [50, 50]
        data.n2 = [100, 100]
    else:
        #Fornax A
        filename = '../Fits_files/FnxA.fits'
        data.x1 = [275, 150]
        data.x2 = [825, 900]
        data.n1 = [225, 100]
        data.n2 = [275, 150]

    fits.filename = filename
    fits.getFitsFile()

    data.dataset = fits.source_data
    data.angScale = [fits.xAngScale, fits.yAngScale]
    data.position = [fits.RA, fits.Dec]
    data.pxl_position = [fits.rowDecpxl, fits.colRApxl]

    if plot_zoomed_image:
        dataset = zoomed_image()
        pretty_plot(dataset)

    if plot_image:
        straight_image()
        pretty_plot(data.dataset)

    if plot_noise:
        noise = noise_image(data.dataset)
        pretty_plot(noise)

    if plot_noise_hist:
        noise = noise_image(data.dataset)
        histogram_plot(noise, 'auto')

    if do_circular_fit:
        dataset = zoomed_image()
        model, resids = circular_fit(dataset)
        pretty_plot(dataset)
        pretty_plot(model)
        pretty_plot(resids)

    if plot_circular_moments:
        dataset = zoomed_image()
        model, resids = circular_fit(dataset)
        moments = shapes.makeMomentArray()
        moment_map(moments)

    if do_model_noise_cf:
        straight_image()
        model, resids = circular_fit(data.dataset)
        noise = noise_image(data.dataset)
        model_noise = noise_image(model)
        histogram_2plot(noise, model_noise, 20)

    if do_elliptical_fit:
        dataset = zoomed_image()
        model, resids = elliptical_fit(dataset)
        pretty_plot(dataset)
        pretty_plot(model)
        pretty_plot(resids)

    if plot_elliptical_moments:
        dataset = zoomed_image()
        model, resids = elliptical_fit(dataset)
        moments = shapes.makeMomentArray()
        moment_map(moments)

    if do_model_noise_cf2:
        straight_image()
        model, resids = elliptical_fit(data.dataset)
        noise = noise_image(data.dataset)
        model_noise = noise_image(model)
        histogram_2plot(noise, model_noise, 20)

    if plot_nmse_circular:
        noise = noise_image(data.dataset)
        dataset = zoomed_image()
        model, resids = circular_fit(dataset)
        data.debug = model
        add_moments(dataset, noise)







