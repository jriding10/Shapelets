% calculates signal statistics: normalised mean squared error, PSNR and the
% structural similarity (SSIM).

function [NMSE, PSNR, SSIM] = stat_calc(y, x)

K = size(x);
pxl = sqrt(K(1));
norm = 0;
mse = 0;
PSNR = 0;
SSIM = 0;
k1 = 0.01;
k2 = 0.03;

k=0;

for i=1:pxl
    for j=1:pxl
        k=k+1;
        z(k) = y(k) - x(k);
        mse = mse + z(k)*z(k);
        norm = norm + y(k)*y(k);
    end
end

NMSE = mse/norm;

G = max(max(y));

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









