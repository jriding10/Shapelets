% For visibility data: use uvlist to obtain data -> Data_intep to make it
% meaningful
% 23/5/13 - trying to reconstruct the PupA_dump (X13 data) from its coeffs.
% No PA or offsets used

clear all;
close all;

% data columns are: u(lambda). v(lambda). bl. vis(real). vis(complex).
% weight
% PupA_dump files
rawdata = load('../Text/PupA_vis_old.txt', '-ascii');
fnew = load('../Text/PupA_dump_coeffs_norot.txt', '-ascii');
shapes = load('../Text/PupA_dump_shapes.txt', '-ascii');
check = load('../Text/PupA_dump_image.txt', '-ascii');
n_max = 25;

% Constants
c = 3e8;                   % if you don't know, quit physics
radians = pi/180;          % conversion: degrees -> radians
short_bl = 30;             % short baseline (lambda)
obs_freq = 179855000;	   % freq Hz
res = 0.078125*radians;    % shapes file in terms of pxls, not ang size 

shapes(1)=shapes(1)*res;
shapes(2)=shapes(2)*res;
% shapes(4)=shapes(4)*res;
% shapes(5)=shapes(5)*res;
shapes(3)=0;
shapes(4)=0;
shapes(5)=0;

% Removes short baselines
short = 0;
uvdata = sort(rawdata, 3);
while rawdata(1,3) <= short_bl
    rawdata(1,:)=[];
    short = short +1;
end

s = size(uvdata);
num_bl = s(1);
max_bl = max(uvdata(:,3));

vis_model = zeros(num_bl,1);        % contains the model visibilities
vis_split = zeros(num_bl,2);        % contains the real/imag vis
data = zeros(num_bl,1);             % actual uv data
weight = zeros(num_bl,1);           % weight matrix
coords = zeros(num_bl,2);           % uv coords
coords(:,1) = uvdata(:,2);
coords(:,2) = uvdata(:,1);

% Be wary as to the definition of weight <-> inttime : this assumes inttime
wnorm = max(uvdata(:,6));

for i=1:num_bl
    data(i) = uvdata(i,4)+1i*uvdata(i,5);
    weight(i) = uvdata(i,6)/wnorm;
end

% visibilty data
type = 1;
vis_model = col_reconstruct(coords, shapes, fout, type);
for i=1:num_bl
    vis_split(i,1) = real(vis_model(i));
    vis_split(i,2) = imag(vis_model(i));
end
% save('PupA_vis_model.txt', 'vis_split', '-ascii');
 
% to compare -> dft
pxl = 128;
pxl_info = [pxl, pxl];

theta_max = 1*radians;
theta_min = (2*theta_max)/pxl;
tot_pxl = pxl*pxl;
im_coords = zeros(tot_pxl, 2);
k=0;

for i=1:pxl
    for j=1:pxl
        k=k+1;
        im_coords(k,1)=sin(-theta_max+(i-1)*theta_min);
        im_coords(k,2)=sin(-theta_max+(j-1)*theta_min);
    end
end

type = 0;
im_model = col_reconstruct(im_coords, shapes, fout, type);

[thing1, thing2, thing3] = vis_res(coords, im_coords, shapes, fout, n_max, weight, data) 

% save('PupA_dump_imcoords.txt', 'im_coords', '-ascii');
% save('PupA_dump_viscoords.txt','coords', '-ascii');
% save('PupA_dump_coldata.txt', 'actual', '-ascii');
% save('PupA_dump_fullmodel.txt', 'vis_split', '-ascii');

% All results are images 
% uv reconstruction from the fhat coefficients
uvmodel = my_idft(im_coords, coords, vis_model);
figure(1)
contour(real(uvmodel));
title('Reconstruction from UV Model');

% uv reconstruction from data points
vis_orig = my_idft(im_coords, coords, data);
figure(2)
im_from_vis = real(vis_orig);
contour(im_from_vis);
title('Image from Vis Data');

[NMSE, PSNR, SSIM] = stat_calc(im_from_vis, real(uvmodel))

% resids = im_from_vis - real(result);
% figure(7)
% contour(resids)
% title('Resids');