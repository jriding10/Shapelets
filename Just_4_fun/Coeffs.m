% Shapelets for image data
clear all;
close all;

% read in data
% shapes = load('../Text/FornaxA_180M_shapes.txt', '-ascii');
% check2 = load('../Text/FornaxA_180M_data.txt','-ascii');
% f_fixed = load('../Text/FornaxA_180M_n25_fixed.txt', '-ascii');
% coords = load('../Text/Fornax_smallcoords.txt', '-ascii');
% shapes = load('../Text/ForA_180M_cirshapes.txt', '-ascii');
% f_fixed = load('../Text/ForA_180M_circoeffs.txt', '-ascii');
shapes = load('../Text/PupA_shapes2.txt', '-ascii');
check2 = load('../Text/PupA_data.txt','-ascii');
f_fixed = load('../Text/PupA_coeffs2.txt', '-ascii');
coords = load('../Text/PupA_coords.txt', '-ascii')*0.02*pi/180;

coeffs = size(f_fixed);

[sorts, ind] = sort(abs(f_fixed(:,3)), 'descend');

for i=1:coeffs(1)
    num = ind(i);
    fout(i,1) = f_fixed(num,1);
    fout(i,2) = f_fixed(num,2);
    fout(i,3) = f_fixed(num,3);
end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imsize = size(check2);
xbit = imsize(1);
ybit = imsize(2);
pxl_info = [xbit, xbit];

k=0;
fleft = fout;
remaining = 1;

type = 0;
while remaining > 0
    k=k+1;
    fsize = size(fleft);
    fhat = fsize(1);
    model = full_reconstruct(coords, shapes, fleft, type, pxl_info);
    residual = check2 - model;
    [NMSE, PSNR, SSIM] = stat_calc(check2, model);
    study(k,1)=coeffs(1)-fhat;
    study(k,2)=NMSE;
    study(k,3)=PSNR;
    study(k,4)=SSIM;
    new_size = fhat-5
    if new_size > 0
        fleft = zeros(new_size,3);
        for i=1:new_size
            fleft(i,:)=fout(i,:);
            remaining = 0;
        end
    end
    if new_size < 0
        remaining = 0;
    end
end

% save('../Text/PupA_coeffstudy.txt', 'study', '-ascii');



