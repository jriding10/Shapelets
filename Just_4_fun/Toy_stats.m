% calculates signal statistics: normalised mean squared error, PSNR and the
% structural similarity (SSIM).

function [NMSE, PSNR, SSIM] = stat_calc(data, model)

K = size(model);
row_pxl = K(1);
col_pxl = K(2);
norm = 0;
mse = 0;
PSNR = 0;
SSIM = 0;
k1 = 0.01;
k2 = 0.03;
n = zeros(row_pxl*col_pxl,1);

k=0;

residual = data - model;

for i=1:row_pxl
    for j=1:col_pxl
        k=k+1;
        mse = mse + residual(i,j)*residual(i,j);
        norm = norm + data(i,j)*data(i,j);
        x(k) = data(i,j);
        y(k) = model(i,j);
        z(k) = residual(i,j);
%         n(k) = noise(k);
    end
end

NMSE = mse/norm;

G = max(max(data));

PSNR = 20*log(G)-10*log(1/(row_pxl*col_pxl)*mse);

c1 = k1*G;
c2 = k2*G;

mux = mean(x);
muy = mean(y);
covxy = mean(x.*y)-mux*muy;
varx = var(x);
vary = var(y);

SSIM = (2*mux*muy+c1)*(2*covxy+c2)/((mux*mux+muy*muy+c1)*(varx+vary+c2));










