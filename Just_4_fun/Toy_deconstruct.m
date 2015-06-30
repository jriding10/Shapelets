% an image based on stored parameters - designed for different
% resolutions and extents in x and y

function f_coeffs = deconstruct(coords, n_max, shapes, Idata)
beta=[shapes(1), shapes(2)];
len = size(coords);
    
PA = shapes(3);
transform = [cos(PA), -sin(PA); sin(PA), cos(PA)];
shifted_coords(:,1) = coords(:,1)-shapes(4)*ones(len(1),1);
shifted_coords(:,2) = coords(:,2)-shapes(5)*ones(len(1),1);
rot_coords = transform*transpose(shifted_coords);
rot_coords = transpose(rot_coords);

M = Toy_basis(rot_coords, n_max, beta);

i=0;
for n1=0:n_max
    for n2=0:n_max
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
    end
end
