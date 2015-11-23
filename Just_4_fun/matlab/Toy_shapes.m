% Script to do everything in pixel coords - designed to laughs
% Shapelets for image data
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical constants
n_max = 25;                        % maximum shapelet order

% Image specific values
% data = imread('Starship.jpg', 'jpg');
data = imread('Serenity.png', 'png');
% resize
for i=250:850
    for j=600:1300
        data3(i-249, j-599,:) = data(i,j,:);
    end
end

data2 = rgb2gray(data3);
data2 = double(data2);
imsize = size(data2);

row_pxls = imsize(1);
col_pxls = imsize(2);
sqpxls = row_pxls*col_pxls;

coords = zeros(sqpxls,2);
coldata = zeros(sqpxls,1);

k=0;
for i=1:row_pxls
    for j=1:col_pxls
        k=k+1;
        coldata(k) = data2(i,j);
        coords(k,1) = i - round(row_pxls/2);
        coords(k,2) = j - round(col_pxls/2);
    end
end




pxl_info = [row_pxls, col_pxls];
[big, small, PA, offsets] = cov_fit(coords, coldata);

% My PA is taken as the anti-clockwise rotation from the x axis: this works
% major = max(row_pxls, col_pxls);
% minor = min(row_pxls, col_pxls);
major=300;
minor=300;
shapes(1) = 0.9*major*power(n_max, -0.52);
shapes(2) = 0.9*minor*power(n_max, -0.52);
shapes(3) = -PA;
shapes(4) = offsets(1);
shapes(5) = offsets(2);
% shapes(3) = 0;
% shapes(4) = 0;
% shapes(5) = 0;

fout = Toy_deconstruct(coords, n_max, shapes, coldata);
model = Toy_reconstruct(coords, shapes, fout, pxl_info);
[NMSE, PSNR, SSIM] = Toy_stats(data2, model)

residual = data2-model;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
surf(data2);
view(0,90);
figure(2)
surf(model);
view(0,90);
figure(3)
surf(residual);
view(0,90);

% 
% out = untangle(fnew);
% % 
% n = 0:n_max;
% figure(4)
% surf(n,n,out);
% view(0,90);

