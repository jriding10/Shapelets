#! /bin/bash python
# Noise.py

##########################################################################
## Does the full shapelet build etc but without the gui interface. The  ##
## idea is to examine the results in terms of noise performance.        ##
##########################################################################

import numpy as np
import matplotlib.pyplot as plt

import Fits
import Source
import Coords
import Shapes


# Pick a source: PKS0410, 3C353, ForA
which_source = 'ForA'
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
plot_nmse_circular = False
plot_nmse_elliptical = False
plot_circular_beta = False
set_circular_beta = False
plot_circular_nmax = False
plot_circular_moments_ssim = False
use_x_moments = False
do_gaussian_fit = False
do_gaussian_series_fit = True
do_egauss_fit = False
do_egauss_series_fit = False

class Common:
    def __init__(self):
        self.nmax = 15
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
    data.position = coords.changeRADecRef(data.pxl_position, data.x1, data.position, data.angScale)
    data.pxl_position = data.x1
    dataset = coords.resizeData(data.dataset, data.x1, data.x2)
    data.coord_list = coords.calcCoordinateList(data.angScale)
    data.raAxis = coords.raCoords(data.position, data.angScale)
    data.decAxis = coords.decCoords(data.position, data.angScale)
    return dataset

def offset_coords(dataset):
    source.extended_source = dataset
    weighted_centre = source.calcWeightedCentre(data.coord_list)
    data.coord_list[:, 0] -= weighted_centre[0]
    data.coord_list[:, 1] -= weighted_centre[1]

def offset_coords_maxflux(dataset):
    source.extended_source = dataset
    centre = source.calcMaxFluxCentre(data.coord_list)
    data.coord_list[:, 0] -= centre[0]
    data.coord_list[:, 1] -= centre[1]

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

def default_fit(dataset):
    extent = find_extent()
    source.extended_source = dataset
    major, minor, pa = source.calcDefaultValues(data.angScale, extent)
    covar_matrix = source.calcFluxMatrix(extent, data.coord_list)
    eigenval, eigenvect = source.calcEigenvalueDecomp(covar_matrix)
    pa = source.calcPositionAngle(eigenvect)
    shapes.position_angle = pa
    coords = source.calcTransformCoordinates(shapes.position_angle, data.coord_list)
    data.coord_list = coords
    model, resids = beta_fit(dataset, major, minor)
    return model, resids

def circular_default_fit(dataset):
    extent = find_extent()
    major, minor, pa = source.calcDefaultValues(data.angScale, extent)
    model, resids = beta_fit(dataset, minor, minor)
    return model, resids

def circular_fit(dataset, beta):
    shapes.major = beta
    shapes.minor = beta
    shapes.position_angle = 0.0
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

def add_moments(dataset):
    flat_data = dataset.flatten()
    total_moments = shapes.moments.shape[0]
    measure = np.zeros((total_moments, 1))
    for i in range(total_moments):
        measure[i], model = only_use_moments(flat_data, i+1)
    return measure

def only_use_moments(flat_data, num_moments):
    moments = shapes.minimiseMoments(num_moments)
    flat_model = shapes.calcShapeletModel(data.coord_list, moments)
    model = coords.expandArray(flat_model)
    new_model = model.flatten()
    ssim = source.calcSSIM(flat_data, new_model)
    return ssim, model

def changing_beta(dataset):
    side_length = max(dataset.shape[0], dataset.shape[1])
    flat_data = dataset.flatten()
    min_beta = min(fits.yAngScale, fits.xAngScale)
    if data.nmax != 0:
        max_beta = 0.9*side_length*min_beta*np.power(data.nmax, -0.52)
    else:
        max_beta = 0.5*side_length*min_beta
    num_points = int((max_beta-min_beta)/min_beta)
    measures = np.zeros((num_points, 3))
    beta = min_beta
    for i in range(num_points):
        major = beta + i*min_beta
        model, resids = circular_fit(dataset, major)
        flat_model = model.flatten()
        measures[i, 0] = np.degrees(shapes.major-min_beta)*60
        measures[i, 1] = source.calcNMSE(flat_data, flat_model)
        measures[i, 2] = source.calcSSIM(flat_data, flat_model)
    return measures

def find_extent():
    extent = [0, 0]
    extent[0] = int(data.x2[0] - data.x1[0])
    extent[1] = int(data.x2[1] - data.x1[0])
    return extent

def beta_fit(dataset, major, minor):
    shapes.major = major
    shapes.minor = minor
    model, resids = create_shapelet_fit(dataset)
    return model, resids

def change_nmax(dataset, major, minor):
    flat_data = dataset.flatten()
    measures = np.zeros((26, 2))
    for i in range(25):
        data.nmax = i
        model, resids = beta_fit(dataset, major, minor)
        flat_model = model.flatten()
        measures[i, 0] = i
        measures[i, 1] = source.calcSSIM(flat_data, flat_model)
    return measures

def calc_gaussian(dataset):
    data.nmax = 0
    measures = changing_beta(dataset)
    maxSSIM = max(measures[:, 2])
    k = 0
    while measures[k, 2] < (maxSSIM-0.00001):
        k += 1

    return measures[k, :]

def calc_egauss(dataset):
    min_beta = data.angScale[0]
    model, resids = default_fit(dataset)
    coords = source.calcTransformCoordinates(shapes.position_angle, data.coord_list)
    data.coord_list = coords
    flat_data = dataset.flatten()
    steps = 20
    data.nmax = 0
    k = 0
    measures = np.zeros((steps*steps, 3))
    for i in range(steps):
        beta1 = min_beta*(i+1)
        for j in range(steps):
            beta2 = min_beta*(j+1)
            model, resids = beta_fit(dataset, beta1, beta2)
            flat_model = model.flatten()
            measures[k, 0] = beta1
            measures[k, 1] = beta2
            measures[k, 2] = source.calcSSIM(flat_data, flat_model)
            k += 1

    k = 0
    best = np.max(measures[:,2])
    while measures[k, 2] < best:
        k += 1
    return measures[k, :]

def gaussian_series(dataset):
    first_gauss = calc_gaussian(dataset)
    flat_data = dataset.flatten()
    good_fit = first_gauss[2]
    beta = np.radians(first_gauss[0]/60)
    model, resids = beta_fit(dataset, beta, beta)
    flat_model = model.flatten()
    good_fit = source.calcSSIM(flat_data, flat_model)
    k = 0

    while good_fit < 0.995:
        k+=1
        print("Gaussian number %f." % k)
        offset_coords_maxflux(resids)
        opt_fit = calc_gaussian(resids)
        beta = np.radians(opt_fit[0]/60)
        if beta <= np.radians(1.0/3600):
            good_fit = 1.0
        else:
            model_piece, resids = beta_fit(resids, beta, beta)
            model += model_piece
            data.debug2 = model
            flat_model = model.flatten()
            good_fit = source.calcSSIM(flat_data, flat_model)
    flat_resids = flat_data - flat_model
    nmse = source.calcNMSE(flat_data, flat_resids)
    resids = coords.expandArray(flat_resids)
    print("There are %f gaussians  with a nmse of %f and a ssim of %f." % (k, nmse, good_fit))
    return model, resids

def egauss_series(dataset):
    first_egauss = calc_egauss(dataset)
    flat_data = dataset.flatten()
    good_fit = first_egauss[2]
    model, resids = beta_fit(dataset, first_egauss[0], first_egauss[1])
    flat_model = model.flatten()
    k = 0

    while good_fit < 0.995:
        k+=1
        offset_coords_maxflux(resids)
        opt_fit = calc_egauss(resids)
        if min(opt_fit[0], opt_fit[1]) <= np.radians(1.0/3600):
            good_fit = 1.0
        else:
            model_piece, resids = beta_fit(resids, opt_fit[0], opt_fit[1])
            model += model_piece
            flat_model = model.flatten()
            good_fit = source.calcSSIM(flat_data, flat_model)
    flat_resids = flat_data - flat_model
    nmse = source.calcNMSE(flat_data, flat_resids)
    resids = coords.expandArray(flat_resids)
    print("There are %f gaussians with a nmse of %f and a ssim of %f." % (k, nmse, good_fit))
    return model, resids

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
        offset_coords(dataset)
        model, resids = circular_default_fit(dataset)
        pretty_plot(dataset)
        pretty_plot(model)
        pretty_plot(resids)

    if plot_circular_moments:
        dataset = zoomed_image()
        offset_coords(dataset)
        model, resids = circular_default_fit(dataset)
        moments = shapes.makeMomentArray()
        moment_map(moments)

    if do_model_noise_cf:
        straight_image()
        offset_coords(data.dataset)
        model, resids = circular_default_fit(data.dataset)
        noise = noise_image(data.dataset)
        model_noise = noise_image(model)
        histogram_2plot(noise, model_noise, 20)

    if do_elliptical_fit:
        dataset = zoomed_image()
        offset_coords(dataset)
        model, resids = default_fit(dataset)
        pretty_plot(dataset)
        pretty_plot(model)
        pretty_plot(resids)

    if plot_elliptical_moments:
        dataset = zoomed_image()
        offset_coords(dataset)
        model, resids = default_fit(dataset)
        moments = shapes.makeMomentArray()
        moment_map(moments)

    if do_model_noise_cf2:
        straight_image()
        offset_coords(data.dataset)
        model, resids = default_fit(data.dataset)
        noise = noise_image(data.dataset)
        model_noise = noise_image(model)
        histogram_2plot(noise, model_noise, 20)

    if plot_nmse_circular:
        noise = noise_image(data.dataset)
        dataset = zoomed_image()
        offset_coords(dataset)
        model, resids = circular_default_fit(dataset)
        add_moments(dataset, noise)

    if plot_circular_beta:
        dataset = zoomed_image()
        offset_coords(dataset)
        measures = changing_beta(dataset)

    if set_circular_beta:
        dataset = zoomed_image()
        offset_coords(dataset)
        if which_source == 'PKS0410':
            beta = np.radians(0.25/60)
        elif which_source == '3C353':
            beta = np.radians(0.92/60)
        else:
            beta = min(fits.xAngScale, fits.yAngScale) + np.radians(325.0/3600)
        model, resids = beta_fit(dataset, beta, beta)
        moments = shapes.makeMomentArray()
        pretty_plot(dataset)
        pretty_plot(model)
        pretty_plot(resids)
        moment_map(moments)

    if plot_circular_nmax:
        dataset = zoomed_image()
        offset_coords(dataset)
        if which_source == 'PKS0410':
            beta = np.radians(0.25/60)
        elif which_source == '3C353':
            beta = np.radians(0.92/60)
        else:
            beta = min(fits.xAngScale, fits.yAngScale) + np.radians(325.0/3600)
        measures = change_nmax(dataset, beta, beta)
        plt.plot(measures[:,0], measures[:,1])
        plt.xlabel('n_{max}')
        plt.ylabel('SSIM')
        plt.show()

    if plot_circular_moments_ssim:
        dataset = zoomed_image()
        offset_coords(dataset)
        if which_source == 'PKS0410':
            beta = np.radians(0.25/60)
        elif which_source == '3C353':
            beta = np.radians(0.92/60)
        else:
            beta = min(fits.xAngScale, fits.yAngScale) + np.radians(325.0/3600)
        model, resids = circular_fit(dataset, beta)
        ssim = add_moments(dataset)
        plt.plot(ssim)
        plt.xlabel('Number of Moments')
        plt.ylabel('SSIM')
        plt.show()

    if use_x_moments:
        dataset = zoomed_image()
        offset_coords(dataset)
        if which_source == 'PKS0410':
            beta = np.radians(0.25/60)
        elif which_source == '3C353':
            beta = np.radians(0.92/60)
        else:
            beta = min(fits.xAngScale, fits.yAngScale) + np.radians(325.0/3600)
        model, resids = circular_fit(dataset, beta)
        flat_data = dataset.flatten()
        ssim, model = only_use_moments(flat_data, 60)
        pretty_plot(model)

    if do_gaussian_fit:
        dataset = zoomed_image()
        offset_coords(dataset)
        opt_fit = calc_gaussian(dataset)
        beta = np.radians(opt_fit[0]/60)
        model, resids = beta_fit(dataset, beta, beta)
        pretty_plot(dataset)
        pretty_plot(model)
        pretty_plot(resids)

    if do_gaussian_series_fit:
        dataset = zoomed_image()
        offset_coords(dataset)
        model, resids = gaussian_series(dataset)
        pretty_plot(dataset)
        pretty_plot(model)
        pretty_plot(resids)

    if do_egauss_fit:
        dataset = zoomed_image()
        offset_coords(dataset)
        major, minor, ssim = calc_egauss(dataset)
        model, resids = beta_fit(dataset, major, minor)
        print('The ssim is %f' % ssim)
        pretty_plot(dataset)
        pretty_plot(model)
        pretty_plot(resids)

    if do_egauss_series_fit:
        dataset = zoomed_image()
        offset_coords(dataset)
        model, resids = egauss_series(dataset)
        pretty_plot(dataset)
        pretty_plot(model)
        pretty_plot(resids)







