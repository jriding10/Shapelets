% computes the basis function M at the approipate visibility points
% hermite_rec(n) generates the co-effs of a hermite polynomial for eg. H2(x)=4x^2-2, 
% hermite_rec(2) gives 4 0 -2 source:
% http://suinotes.wordpress.com/2010/05/26/hermite-polynomials-with-matlab/


function [M,nvals] = basis_functions(coords, n_max, beta, domain, unique)

n1 = 0;                                  
n2 = 0;
q = 1;				 % index for M(n1,n2)

% domain = 0 is a visibility basis, 1 is an image basis
if domain == 0
	scale = beta;
	active = 1;
else
	scale = 1/(beta);
	active = 0;
end

for n1 = 0:n_max
    h1 = hermite_rec(n1);
    n2 = 0;
    while (n1+n2) <= n_max
        h2 = hermite_rec(n2);
        FTfactor = power(i, n1+n2);
        FTfactor = power(FTfactor, active);
        D = sqrt(power(2,n1+n2)*pi*1/(scale*scale)*factorial(n1)*factorial(n2));

    	for j=1:unique
            baseline = coords(j,3); 
            C(j,q) = FTfactor*exp(-0.5*baseline*baseline*scale*scale);
            k=1;
            h1_ans=0;
            h2_ans=0;
	    % calculates the polynominal (co-effs in h1) at x=u
            while k <= (n1+1)		
            	h1_ans = h1_ans+h1(k)*power(scale*coords(j,1), (n1+1)-k);
            	k=k+1;
            end
            k=1;
	    % calculates the polynomial (co-effs in h2) at x=v
            while k <= (n2+1)		
             	h2_ans = h2_ans+h2(k)*power(scale*coords(j,2), (n2+1)-k);
             	k=k+1;
            end
            Vc_bar(j,q) = C(j,q)/D*h1_ans*h2_ans;
        end

	M(:,n1+1,n2+1) = Vc_bar(:,q);
    nvals(q,1)=n1;
    nvals(q,2)=n2;
	n2 = n2+1;
	q = q+1;
     end
end

