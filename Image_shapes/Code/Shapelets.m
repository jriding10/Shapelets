% Shapelets: Takes a fits files input and outputs the information needed
% for shapelets - betas, moments, PA etc. May include multicomponent fits
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user inputs

fits_file_in = '~/Data/CenA.fits';
existing_beta = 0;      % 1 = using an existing file (ie not generating a new one)
para_file = '../Text/J0351_paras.txt';
existing_fit = 0;       % 1 = using an existing file 
coeff_file = '../Text/J0351_351coeffs.txt';
obj_size_max = 25;      % major axis, arcmins (from NED)

n_max = 25;            % use all the available coefficients = 25
tot_coeffs = 351;        % number of moments to include in fit
n_approx = 5;

% minimise coefficients
min_coeffs = 0;

% save files
save_beta = 1;          % save a parameter file
new_para_file = '../Text/CenA_paras.txt';
save_coeffs = 1;        % save the coefficient file
new_coeff_file = '../Text/CenA_351coeffs.txt';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical constants
radians = pi/180;                   % degrees to radians
n_approx = 5;                       % autobeta speedup
coeffs = 0.5*(n_max+1)*(n_max+2);   % number of moments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading and initialising data and variables

data = fitsread(fits_file_in, 'primary');
info = fitsinfo(fits_file_in);

if existing_beta == 1
    shapes = load(para_file, '-ascii');    
end

if existing_fit == 1
    moments = load(coeff_file, '-ascii');
end

[pix_obj, coord_obj, ang_res]=fits_info(info); 
obj = [coord_obj, obj_size_max];
[precoords, predataset, col_data, ra_and_dec] = prepare_data(pix_obj, obj, ang_res, data);

% Code parameters - image rows = y axis = dec = axis1
data_size = size(predataset);
npix_side = data_size(1);

npix_mid = (npix_side/2);
npix = npix_side*npix_side;

% Prepare data and axes    
if existing_beta == 0  
    [big, small, PA, offsets] = cov_fit(precoords, predataset);
    PA = -PA;
    shapes(1) = (coord_obj(1)*radians + offsets(2))/radians;
    shapes(2) = (coord_obj(2)*radians + offsets(1))/radians;
    shapes(5) = PA/radians;     
end

dec_cen = pix_obj(2)+round((coord_obj(2)-shapes(2))*radians/ang_res);
ra_cen = pix_obj(1)+round((coord_obj(1)-shapes(1))*radians/ang_res);
midpt = [dec_cen, ra_cen];

[coords, dataset, col_data, ra_and_dec] = prepare_data(midpt, obj, ang_res, data);

ra_axis = ra_and_dec(:,1);
dec_axis = ra_and_dec(:,2);
coord_size = size(coords);
npix = coord_size(1);
npix_side = sqrt(npix);

[rough_nmse, check] = estinoise(dataset);

PA = shapes(5)*radians;
transform = [cos(PA), -sin(PA); sin(PA), cos(PA)];
rot_coords = transform*transpose(coords);
im_coords = transpose(rot_coords);

col_mod = zeros(npix,1);		% model in a column
col_resid = zeros(npix,1);		% residuals in a column
model = zeros(sqrt(npix));              % model image
residual = zeros(sqrt(npix));		% residual image

if existing_beta == 0
    res = [ang_res, obj_size_max];
    [beta1, beta2, temp] = majorminor(dataset, n_approx, im_coords, res);
    shapes(3) = beta1*60/radians;
    shapes(4) = beta2*60/radians;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if existing_fit == 0
    moments = deconstruct(im_coords, ang_res, n_max, shapes, col_data);
end

[sorts, ind] = sort(abs(moments(:,3)), 'descend');
for i=1:tot_coeffs
    k = ind(i);
    sorted(k,:) = moments(k,:);
end

if min_coeffs == 1
    [tot_coeffs, study] = minco_long(im_coords, sorted, col_data, shapes, rough_nmse);
end 

for i=1:tot_coeffs
    final_moments(i,:) = sorted(i,:);
end

if save_beta == 1
     save(new_para_file, 'shapes', '-ascii', '-double');
end

if save_coeffs == 1
    save(new_coeff_file, 'final_moments', '-ascii', '-double');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build model and assess performance
if tot_coeffs > 0
    col_mod = reconstruct(im_coords, ang_res, shapes, final_moments);
end

%rescale = inv(transpose(col_mod)*col_mod)*transpose(col_mod)*(col_data);
rescale = 1;

[NMSE, PSNR, SSIM] = stat_calc(col_data, col_mod)
col_resid = col_data - col_mod;

residual = convert(col_resid, npix_side, npix_side);
model = convert(col_mod, npix_side, npix_side);

% normalising a colourbar across images based on data
data_min = min(min(dataset));
data_max = max(max(dataset));

figure(1)
%surf(ra_axis, dec_axis, dataset, 'LineStyle', 'None'); 
%axis image;
%caxis([data_min, data_max]);
%colorbar;
%view(0,90);
contour(ra_axis, dec_axis, dataset);

figure(2)
%surf(ra_axis, dec_axis, model, 'LineStyle', 'None'); 
%axis image;
%caxis([data_min, data_max]);
%colorbar;
%view(0,90);
contour(ra_axis, dec_axis, model);

figure(3)
%surf(ra_axis, dec_axis, residual, 'LineStyle', 'None'); 
%axis image;
%caxis([data_min, data_max]);
%colorbar;
%view(0,90);
contour(ra_axis, dec_axis, residual);

