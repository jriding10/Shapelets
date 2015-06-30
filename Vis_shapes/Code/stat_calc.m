% calculates signal statistics: normalised mean squared error, PSNR and the
% structural similarity (SSIM).

function [NMSE, PSNR, SSIM] = stat_calc(image, model)

K = size(model);
pxl = K(1);
norm = 0;
mse = 0;
PSNR = 0;
SSIM = 0;
k1 = 0.01;
k2 = 0.03;
n = zeros(pxl*pxl,1);

k=0;

residual = image - model;

for i=1:pxl
    for j=1:pxl
        k=k+1;
        mse = mse + residual(i,j)*residual(i,j);
        norm = norm + image(i,j)*image(i,j);
        x(k) = image(i,j);
        y(k) = model(i,j);
        z(k) = residual(i,j);
%         n(k) = noise(k);
    end
end

NMSE = mse/norm;

G = max(max(image));

PSNR = 20*log(G)-10*log(1/(pxl^2)*mse);

c1 = k1*G;
c2 = k2*G;

mux = mean(x);
muy = mean(y);
covxy = mean(x.*y)-mux*muy;
varx = var(x);
vary = var(y);

SSIM = (2*mux*muy+c1)*(2*covxy+c2)/((mux*mux+muy*muy+c1)*(varx+vary+c2));

% mux = mean(z);
% muy = mean(n);
% covxy = mean(z.*n)-mux*muy;
% varx = var(z);
% vary = var(n);
% 
% SSIM2 = (2*mux*muy+c1)*(2*covxy+c2)/((mux*mux+muy*muy+c1)*(varx+vary+c2));









