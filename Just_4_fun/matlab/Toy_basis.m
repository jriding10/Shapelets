% computes the basis functions for shapelets

function [M] = basis(coords, n_max, beta)
hermites_n = load('../Common/hermite_n46_coeffs.txt', '-ascii');
n1 = 0;                                  
n2 = 0;
scale_x=beta(1);
scale_y=beta(2);
csize = size(coords);
tot_coords = csize(1);

for n1 = 0:n_max
    for n2 = 0:n_max
        D = sqrt(power(2,n1+n2)*pi*scale_x*scale_y*factorial(n1)*factorial(n2));

    	for j=1:tot_coords
            C(j) = exp(-0.5*(coords(j,1).^2/(scale_x*scale_x)+coords(j,2).^2/(scale_y*scale_y)));
            k=1;
            h1_ans=0;
            h2_ans=0;
	    % calculates the polynominal (co-effs in h1) at x=u or l
            while k <= (n1+1)		
            	h1_ans = h1_ans+hermites_n(n1+1, k)*power(coords(j,1)/scale_x, (n1+1)-k);
            	k=k+1;
            end
            k=1;
	    % calculates the polynomial (co-effs in h2) at x=v or m
            while k <= (n2+1)		
             	h2_ans = h2_ans+hermites_n(n2+1, k)*power(coords(j,2)/scale_y, (n2+1)-k);
             	k=k+1;
            end
            Vc_bar(j) = C(j)/D*h1_ans*h2_ans;
        end

	M(:,n1+1,n2+1) = Vc_bar(:);
    end
end

