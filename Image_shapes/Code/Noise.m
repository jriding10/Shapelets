% estimates the noise of an image

% Setting the scene...
% physical constants
clear all;
close all;
radians = pi/180;                   % degrees to radians
fraction = 0.1;
fact = 20;                          % scales data to appro. size
raw_data = 0;
numbins = 100;
num_coeffs = 351;
sspan = 1;

% Input file
data_file = '~/Data/ForA.fits';
noise_file = '~/Data/HydA/HydA_noise.txt';
data = fitsread(data_file, 'primary');
info = fitsinfo(data_file);
mom = load('../Text/ForA_351coeffs.txt', '-ascii');
shapes = load('../Text/ForA_paras.txt', '-ascii');
obj_size = 1.0;

[pix_obj, coord_obj, ang_res]=fits_info(info); 

if num_coeffs > 0
    for i=1:num_coeffs
        moments(i,:) = mom(i,:);
    end
end

obj = [coord_obj, obj_size];
dec_cen = pix_obj(2)+round((coord_obj(2)-shapes(2))*radians/ang_res);
ra_cen = pix_obj(1)+round((coord_obj(1)-shapes(1))*radians/ang_res);
midpix = [dec_cen, ra_cen];
[coords, dataset, col_data, skycoords] = prepare_data(midpix, obj, ang_res, data);

coord_size = size(coords);
npix = coord_size(1);
npix_side = sqrt(npix);
data_max = max(max(data));

% Initilise matrices
im_coords = zeros(npix,2);            % creates coords for image          
col_mod = zeros(npix,1);
model = zeros(npix_side);             % fit

% Prepare data and axes

PA = shapes(5)*radians;
transform = [cos(PA), -sin(PA); sin(PA), cos(PA)];
rot_coords = transform*transpose(coords);
im_coords = transpose(rot_coords);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise estimate near the source

delta = npix_side*fraction;
ntot = 0;
k=0;
for i=1:npix_side
    for j=1:npix_side
       if i < delta
           k=k+1;
           noise(k)=data(i,j);
           cutout(i,j)=data(i,j);
           ntot = ntot+noise(k);
        end
        if i > (npix_side-delta)
            k=k+1;
            noise(k)=data(i,j);
            cutout(i,j)=data(i,j);
            ntot = ntot+noise(k);
        end
        if j < delta && i > delta && i < (npix_side - delta)
            k=k+1;
            noise(k)=data(i,j);
            cutout(i,j)=data(i,j);
            ntot = ntot+noise(k);
        end
        if j > (npix_side-delta) && i > delta && i < (npix_side - delta)
            k=k+1;
            noise(k)=data(i,j);
            cutout(i,j)=data(i,j);
            ntot = ntot+noise(k);
        end
    end
end

if num_coeffs > 0
    col_mod = reconstruct(im_coords, ang_res, shapes, moments);
end
col_resid = col_data - col_mod;

%% Output stats
min_noise = min(noise);
max_noise = max(noise);
min_resid = min(col_resid);
max_resid = max(col_resid);

noise_avg = mean(noise);
noise_vari = var(noise);
resid_avg = mean(col_resid);
resid_vari = var(col_resid);

nnpix = size(noise);

[NMSE, PSNR, SSIM] = stat_calc(col_data, col_mod)
% generates gaussian noise
cf_noise = noise_avg + sqrt(noise_vari)*randn(nnpix(2),1);

% generate binning vector
% min_bin = min(min_noise, min_resid);
% max_bin = max(max_noise, max_resid);
% bin_size = (max_bin - min_bin)/numbins;
% bin_vec(1) = min_bin;
bin_size = (max_noise-min_noise)/numbins;
bin_vec(1) = min_noise;

for i=2:numbins
    bin_vec(i) = bin_vec(i-1)+bin_size;
end

for i=1:(numbins)
    x = bin_vec(i);
    nnorm = sqrt(2*pi*noise_vari);
    ngauss = exp(-1*(x*x/(2*noise_vari)));
    white(i) = 1/nnorm*ngauss;
    rnorm =  sqrt(2*pi*resid_vari);
    rgauss = exp(-1*(x*x/(2*resid_vari)));
    clean(i) = 1/rnorm*rgauss;
end

hist_noise = hist(noise, bin_vec);
hist_resids = hist(col_resid, bin_vec);

hist_noise_smooth = smooth(hist_noise, sspan);
hist_resids_smooth = smooth(hist_resids, sspan);
 
% figure(1)
% hold on;
% plot(bin_vec, hist_noise_smooth);
% plot(bin_vec, hist_resids_smooth, 'r');
% legend('Intrinsic', 'Residual');
% hold off;
 
for i=1:(numbins-2*sspan)
    new_vec(i) = bin_vec(i+sspan);
    new_noise(i) = hist_noise_smooth(i+sspan);
    new_resids(i) = hist_resids_smooth(i+sspan);
    gnoise(i) = white(i+sspan);
    gresids(i) = clean(i+sspan);
end

% figure(2)
% hold on;
% plot(new_vec, new_noise);
% plot(new_vec, new_resids, 'r');
% legend('Intrinsic', 'Residual');
% hold off;

nscale = inv(gnoise*transpose(gnoise))*(gnoise)*transpose(new_noise);
rscale = inv(gresids*transpose(gresids))*(gresids)*transpose(new_resids);

% figure(3)
% hold on;
% plot(new_vec, new_noise, 'r');
% plot(new_vec, nscale*gnoise);
% legend('Noise', 'Fitted WGN');
% hold off;

figure(2)
hold on;
plot(new_vec, new_resids, 'r');
plot(new_vec, rscale*gresids);
plot(new_vec, nscale*gnoise, 'k');
legend('Residuals', 'Fitted WGN');
hold off;

rout = 0;
for i=1:sspan
    rout = rout+hist_resids(i);
    rout = rout+hist_resids(i+numbins-sspan);
end

rinner = 0;
for i=sspan:(numbins-2*sspan+1)
    rinner = rinner+hist_resids(i);
end

nout = 0;
for i=1:sspan
    nout = nout+hist_noise(i);
    nout = nout+hist_noise(i+numbins-sspan);
end

ninner = 0;
for i=sspan:(numbins-2*sspan+1)
    ninner = ninner+hist_noise(i);
end


% How well does gresid fit the actual data
rMSE = 0;
nMSE = 0;

for i=1:(numbins-2*sspan)
    nMSE = nMSE + (new_noise(i) - nscale*gnoise(i))*(new_noise(i) - nscale*gnoise(i));
    rMSE = rMSE + (new_resids(i) - rscale*gresids(i))*(new_resids(i) - rscale*gresids(i));
end


