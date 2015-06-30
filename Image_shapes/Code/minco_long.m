% finds the point when residuales match the noise

function [needed, temp] = minico(coords, coeffs, col_data, shapes, nmse)

%nmse=rough_nmse;
%coeffs=sorted;
%coords=im_coords;

radians = pi/180;
data_size = size(col_data);
npix = data_size(1);
beta1 = shapes(3)/60*radians;
beta2 = shapes(4)/60*radians;
beta = [beta1, beta2];
col_mod = zeros(npix,1);
col_resid = zeros(npix,1);
factor = 2;

none = max(coeffs(:,1));
ntwo = max(coeffs(:,2));
n_max = max(none, ntwo);
max_flux = max(col_data);

% calculates all the basis functions at once
M = basis(coords, n_max, beta);

resid_nmse = 100;
nmse=nmse/factor;
psnr = 0;
i=0;

% nmse/10 for a margin
%while (nmse<resid_nmse & i<351)
while (psnr<120 & i<351)
    i=i+1;
    n1 = coeffs(i,1);
    n2 = coeffs(i,2);
    fhat = coeffs(i,3);

    Mpiece = M(:,n1+1,n2+1);

    col_mod = col_mod + fhat*Mpiece;
    col_resid = col_data - col_mod;
    mse = 1/npix*sum(col_resid.^2);
    psnr = 20*log(max_flux)-10*log(mse);
    resid_nmse = var(col_resid);
    temp(i,2)=psnr;
    temp(i,1)=resid_nmse;
end

needed = i;
