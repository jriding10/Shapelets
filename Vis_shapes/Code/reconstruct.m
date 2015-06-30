% Creates the (u,v) data model for the source based on pre-computed
% coefficients and a generic basis file

function [vis_data] = reconstruct(baselines, res, B, f_coeffs);
basis = load('basis_fn.txt', '-ascii');

% Program parameters - matrix sizes
bls = size(baselines);
num_bls = bls(1);
bsize = size(basis);
max_size = bsize(2);
fc = size(f_coeffs);
num_coeffs = fc(1);

beta1 = B(1);
beta2 = B(2);
PA = B(3);

% (x,y) are the rotated coordinates of (u_s, v_s), the baseline coords
% stored coeffs take the form n1, n2, coeff

for i=1:num_coeffs
    Mpiece1 = basis(f_coeffs(i,1), :);
    Mpiece2 = basis(f_coeffs(i,2), :);
    f_hat = f_coeffs(i,3);
    for i=1:num_bls
        x1 = (x*beta1)/(base_res*sqrt(2*pi))+1;
        index1 = floor(x1);
        y1low = Mpiece1(index1);                
        y1high = Mpiece1(index1+1);                
        m1 = (y1high - y1low);
        c1 = y1low-m1*index1;  
        thi1=m1*x1+c1;               
        x2 = (rot_coords(i,2)/beta2+mid)/base_res+1;
        index2 = floor(x2);
        y2low = Mpiece2(index2);
        y2high = Mpiece2(index2+1);
        m2 = (y2high - y2low);
        c2 = y2low-m2*index2;
        thi2=m2*x2+c2;
        new_base(i) = 1/sqrt(beta1)*1/sqrt(beta2)*thi1*thi2;
        new_im(i) = new_im(i)+fhat*new_base(i);
     end
end


k=0;
for i=1:pts
    for j=1:pts
        k=k+1;
        final_im(i,j) = new_im(k);
    end
end