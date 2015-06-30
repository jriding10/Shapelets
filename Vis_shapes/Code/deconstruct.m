% an image based on stored parameters - designed for different
% resolutions and extents in x and y

function f_coeffs = deconstruct(coords, n_max, shapes, Idata, type)
% clear all;
% close all;
% coords = load('../Text/FornaxA_coords.txt', '-ascii');
% shapes = load('../Text/FornaxA_shapes.txt', '-ascii');
% Idata = load('../Text/FornaxA_coldata_150M.txt', '-ascii');

len = size(coords);
midpt = sqrt(len(1))/2;
radians = pi/180;
pxl = sqrt(len(1));

beta=[shapes(1), shapes(2)];
    
PA = shapes(3);
transform = [cos(PA), -sin(PA); sin(PA), cos(PA)];
shifted_coords(:,1) = coords(:,1)-shapes(4)*ones(len(1),1);
shifted_coords(:,2) = coords(:,2)-shapes(5)*ones(len(1),1);
rot_coords = transform*transpose(shifted_coords);
rot_coords = transpose(rot_coords);

M = basis(rot_coords, n_max, beta, type);

i=0;
for n1=0:n_max
    n2=0;
    while (n1+n2)<=n_max
        Mpiece = M(:,n1+1,n2+1);
        demon = transpose(Mpiece)*(Mpiece);
        numer = transpose(Mpiece)*Idata;
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

% M00 = convert(M(:,1,1), pxl, pxl);
% contour(f_coeffs(1,3)*M00);