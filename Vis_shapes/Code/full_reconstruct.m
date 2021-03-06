% creates the basis functions, does not use the generic basis file

% function model = full_reconstruct(coords, shapes, f_coeffs, type)
clear all;
close all;

coords = load('PupA_dump_viscoords.txt' , '-ascii');
shapes = load('PupA_dump_shapes.txt', '-ascii');
f_coeffs = load('PupA_dump_coeffs.txt', '-ascii');
imcoords = load('PupA_dump_imcoords.txt', '-ascii');
type = 1;

len = size(coords);
co = size(f_coeffs);
none = max(f_coeffs(:,1));
ntwo = max(f_coeffs(:,2));
n_max = max(none, ntwo);

% Rotation in the image is calculated the same way in the fourier domain,
% so coords can be either (u,v) or (l,m)
PA = shapes(3);
transform = [cos(PA), -sin(PA); sin(PA), cos(PA)];
shifted_coords(:,1) = coords(:,1)-shapes(4)*ones(len(1),1);
shifted_coords(:,2) = coords(:,2)-shapes(5)*ones(len(1),1);
rot_coords = transform*transpose(shifted_coords);
coords = transpose(rot_coords);

beta_x = shapes(1);
beta_y = shapes(2);

M = basis(coords, beta_x, beta_y, type, n_max);

model = zeros(len(1),1);

for i=1:co(1)
    n1 = f_coeffs(i,1);
    n2 = f_coeffs(i,2);
    f_hat = f_coeffs(i,3);
    Mpiece = M(:,n1+1,n2+1);
    model = model+f_hat*(Mpiece);
end

for i=1:len(1)
    if abs(model(i)) < 10e-15
	model(i) = 0.0;  
    end
end

out = my_dft(imcoords, coords, model);