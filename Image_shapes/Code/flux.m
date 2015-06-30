% estimates the contribution of flux from each coefficent

% Setting the scene...
% physical constants
clear all;
close all;
radians = pi/180;                   % degrees to radians
delta = 20;
fact = 50;                          % scales data to appro. size

% Input file
moments = load('../Text/HydraA_351coeffs.txt', '-ascii');
shapes = load('../Text/HydraA_paras.txt', '-ascii');

beta1 = shapes(3)/radians/60;
beta2 = shapes(4)/radians/60;

none = max(moments(:,1));
ntwo = max(moments(:,2));
nmax = max(none, ntwo);
num_coeffs = 0.5*(none+1)*(none+2);

rel_coeffs = floor(nmax/2)+1;
int_flux = zeros(rel_coeffs,2);
rel_moms = zeros(rel_coeffs,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0;
for i=1:num_coeffs
    n1 = moments(i,1);
    n2 = moments(i,2);
    if n1 == 0 || n2 == 0
        if n1==0 && n2==0
            k=k+1;
            rel_moms(k,:)=moments(i,:);
        elseif n1 ~= 0 && mod(n1,2)==0
            k=k+1;
            rel_moms(k,:)=moments(i,:);
        elseif n2 ~= 0 && mod(n2,2)==0
            k=k+1;
            rel_moms(k,:)=moments(i,:);
        end
           
    elseif mod(n1,2) == 0 && mod(n2,2) == 0
        k=k+1;
        rel_moms(k,:) = moments(i,:);
    end
end

factor = sqrt(pi)*sqrt(beta1)*sqrt(beta2);

for i=1:rel_coeffs
    n1 = rel_moms(i,1);
    n2 = rel_moms(i,2);
    nmax = max(n1, n2);
    coeffs = 0.5*(nmax+1)*(nmax+2);
    f12 = rel_moms(i,3);
    upstairs = 1/2*(2-n1-n2);
    c1 = nchoosek(n1, n1/2);
    c2 = nchoosek(n2, n2/2);
    int_flux(i+1,1) = coeffs;
    int_flux(i+1,2) = int_flux(i)+power(2, upstairs)*sqrt(c1)*sqrt(c2)*f12;
end

int_flux = int_flux*factor;

