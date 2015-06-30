% Fits moments to an image and then outputs a fit. Unlike
% deconstruct/reconstruct, recycles basis functions

function NMSE = shape_models(coords, n_max, beta1, beta2, col_data)
%coords = im_coords;
%n_max = 5;
npix = size(coords);
pxl_side = sqrt(npix(1));
radians = pi/180;

beta = [beta1, beta2];
M = basis(coords, n_max, beta);
norm = transpose(col_data)*col_data;

i=0;
for n1=0:n_max
    n2=0;
    while (n1+n2)<=n_max
        Mpiece = M(:,n1+1,n2+1);
        demon = transpose(Mpiece)*(Mpiece);
        numer = transpose(Mpiece)*col_data;
        f = numer/(demon);
        if (n1+n2) == 0
            fmax = f;
        end
        if abs(f) > 0*fmax
            i=i+1;
            f_coeffs(i,1)=n1;
            f_coeffs(i,2)=n2;
            f_coeffs(i,3)=f;
        end
        n2=n2+1;
    end
end

col_model = zeros(npix(1),1);
co = size(f_coeffs);

for i=1:co(1)
    n1 = f_coeffs(i,1);
    n2 = f_coeffs(i,2);
    f_hat = f_coeffs(i,3);
    Mpiece = M(:,n1+1,n2+1);
    col_model = col_model+f_hat*(Mpiece);
end

for i=1:npix(1)
    if abs(col_model(i)) < 10e-15
	col_model(i) = 0.0;  
    end
end

NMSE = 0;
for i=1:npix(1)
    NMSE = NMSE + (col_data(i) - col_model(i))*(col_data(i) - col_model(i));
end
NMSE = 1/norm*NMSE;
