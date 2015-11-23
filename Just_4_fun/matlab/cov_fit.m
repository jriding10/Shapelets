% Weighted covariance matrix to find major, minor and PA
function [major, minor, PA, offsets] = cov_fit(coords, Idata);
% Idata = load('../Text/FornaxA_B00.txt','-ascii');
% coords = load('../Text/FornaxA_coords.txt','-ascii');
% offsets = load('Test_offsets.txt', '-ascii');

S = sum(sum(Idata));
len = size(Idata);
m = len(1)*len(2);
offsets = [0.0,0.0];

x=0;
y = 0;
Sxx = 0;
Sxy = 0;
Syy = 0;

k=0;
for i=1:len(1)
   for j=1:len(2)
	k=k+1;
	x = x+Idata(i,j)*coords(k,1);
	y = y+Idata(i,j)*coords(k,2);
   end
end

x0 = round(x/S);
y0 = round(y/S);

offsets = [x0, y0];

coords(:,1)=coords(:,1)-x0*ones(m,1);
coords(:,2)=coords(:,2)-y0*ones(m,1);

k=0;
for i=1:len(1)
    for j=1:len(2)
        k=k+1;
        Sxx = Sxx + Idata(i,j)*coords(k,1)*coords(k,1);
        Sxy = Sxy + Idata(i,j)*coords(k,1)*coords(k,2);
        Syy = Syy + Idata(i,j)*coords(k,2)*coords(k,2);
    end
end

a11 = Sxx/S;
a12 = Sxy/S;
a22 = Syy/S;

C = [a11, a12; a12, a22];

% Finding the eigenvalue decomposition
[V, lambda] = eig(C);

minor = sqrt(lambda(1,1));
major = sqrt(lambda(2,2));

PA = atan(V(2,2)/V(1,2));




