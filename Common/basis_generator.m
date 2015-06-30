% computes the basis functions for shapelets
clear all;
close all;

hermites_n = load('hermite_coeffs.txt', '-ascii');
n_max = 30;
beta = 1;
res = 0.01;
FoV = 100;
pxl_max = floor(FoV/res)+1;
M = zeros(n_max, pxl_max);
B = zeros(1,pxl_max);
gauss = zeros(pxl_max,1);
coords = zeros(pxl_max,1);
M = zeros(n_max, pxl_max);
B = zeros(1,pxl_max);
gauss = zeros(pxl_max,1);
coords = zeros(pxl_max,1);
radians = pi/180;
h1_ans = zeros(n_max+1,pxl_max);

for i=1:pxl_max
    coords(i)=i*res - FoV/2;
    gauss(i) = exp(-0.5*(coords(i,1)/beta)^2);
end

for n = 0:n_max
        norm = sqrt(power(2,n)*factorial(n)*beta*pi);
        k=1;
    	for j=1:pxl_max
            while k <= (n+1)
            	h1_ans(n+1,j) = h1_ans(n+1,j)+hermites_n(n+1, k)*power(coords(j,1)/beta, (n+1)-k);
                k=k+1;
            end
            k=1;
            B(j) = gauss(j)/norm*h1_ans(n+1,j);
        end

    M(n+1, :) = B(:);
end

save('basis_fn.dat', 'M', '-ascii', '-double');

