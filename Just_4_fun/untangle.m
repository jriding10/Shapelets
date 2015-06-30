% creates a matrix of coefficients

function something = untangle(coeffs)

n1 = max(coeffs(:,1));
n2 = max(coeffs(:,2));
n_max = max(n1,n2);
len = size(coeffs);

something = zeros(n_max+1);

for i=1:len(1)
    n1 = coeffs(i,1);
    n2 = coeffs(i,2);
    something(n1+1, n2+1) = coeffs(i,3);
    relevent(n1+1, n2+1) = 1; 
end

figure(11)
surf(something);
view(0,90);

figure(12)
surf(relevent);
view(0,90);


