% creates the basis functions, does not use the generic basis file

function model = reconstruct(coords, ang_res, shapes, f_coeffs);
% clear all;
% close all;
% f_coeffs = load('../Text/HydraA_351coeffs.txt', '-ascii');
% shapes = load('../Text/HydraA_paras.txt', '-ascii');
% npix = 2809;
% ang_res = 1/60*pi/180;

coeffs = size(f_coeffs);
none = max(f_coeffs(:,1));
ntwo = max(f_coeffs(:,2));
n_max = max(none, ntwo);
npix = size(coords);
nrows = sqrt(npix(1));
midpt = round(nrows/2);
radians = pi/180;

model = zeros(npix(1),1);

beta=[shapes(3), shapes(4)]*radians/60;

M = basis(coords, n_max, beta);

for i=1:coeffs(1)
    n1 = f_coeffs(i,1);
    n2 = f_coeffs(i,2);
    f_hat = f_coeffs(i,3);
    Mpiece = M(:,n1+1,n2+1);
    model = model+f_hat*(Mpiece);
end

for i=1:nrows
    if abs(model(i)) < 10e-15
	model(i) = 0.0;  
    end
end

