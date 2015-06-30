% computes the basis functions for shapelets

function [M] = basis(coords, shapes, type, n_max)
% clear all;
% close all;
% 
% coords = load('../Text/Vis_coords_tester.txt', '-ascii');
% shapes = load('../Text/PupA_dump_cirshapes.txt', '-ascii');
% 
% type=1;
% n_max=25;
% len=size(coords);
% PA = shapes(3);
% n_max=0;
% shapes = pi/180*0.078125*shapes;

% transform = [cos(PA), -sin(PA); sin(PA), cos(PA)];
% shifted_coords(:,1) = coords(:,1)-shapes(4)*ones(len(1),1);
% shifted_coords(:,2) = coords(:,2)-shapes(5)*ones(len(1),1);
% rot_coords = transform*transpose(shifted_coords);
% coords = transpose(rot_coords);

hermites_n = load('../Text/hermite_coeffs.txt', '-ascii');
n1 = 0;                                  
n2 = 0;
csize = size(coords);
tot_coords = csize(1);
beta_x = shapes(1);
beta_y = shapes(2);

% % type 0 is image domain, 1 is visibility 
if type == 0
    shiftl = 0;
    shiftm = 0;
    scale_x = beta_x;
    scale_y = beta_y;     
    xoff = shapes(4);
    yoff = shapes(5);      
end

if type == 1
    scale_x = 1/(2*pi*beta_x);
    scale_y = 1/(2*pi*beta_y);
    shiftl = shapes(4);
    shiftm = shapes(5);
    xoff = 0;
    yoff = 0;
end

for n1 = 0:n_max
    n2 = 0;
    while (n1+n2) <= n_max
        FTfactor = power(1i,n1+n2)^type;
        D = sqrt(power(2,n1+n2)*pi*scale_x*scale_y*factorial(n1)*factorial(n2));

    	for j=1:tot_coords
            x = coords(j,1);
            y = coords(j,2);
            C(j) = exp(-0.5*(x*x/(scale_x*scale_x)+y*y/(scale_y*scale_y)));
            k=1;
            h1_ans=0;
            h2_ans=0;
	    % calculates the polynominal (co-effs in h1) at x=u or l
            while k <= (n1+1)		
            	h1_ans = h1_ans+hermites_n(n1+1, k)*power(x/scale_x, (n1+1)-k);
            	k=k+1;
            end
            k=1;
	    % calculates the polynomial (co-effs in h2) at x=v or m
            while k <= (n2+1)		
             	h2_ans = h2_ans+hermites_n(n2+1, k)*power(y/scale_y, (n2+1)-k);
             	k=k+1;
            end
            shift_term = exp(-1i*2*pi*(shiftl*coords(j,1)+shiftm*coords(j,2)));
%             shift_term = 1;
            Vc_bar(j) = shift_term*FTfactor*C(j)/D*h1_ans*h2_ans; 
%             check(x+uv_start,y+uv_start, n1+1, n2+1) = Vc_bar(j);
        end

	M(:,n1+1,n2+1) = Vc_bar(:);
	n2 = n2+1;
    end
end
    

