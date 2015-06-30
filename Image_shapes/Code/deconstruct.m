% an image based on stored parameters - designed for different
% resolutions and extents in x and y

function f_coeffs = deconstruct(coords, ang_res, n_max, shapes, col_data)

coords_size = size(coords);
npix = coords_size(1);
nrows = sqrt(npix);
midpt = round(nrows/2);
radians = pi/180;

beta=[shapes(3), shapes(4)]*radians/60;

M = basis(coords, n_max, beta);

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
