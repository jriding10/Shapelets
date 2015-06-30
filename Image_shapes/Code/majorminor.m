% Attempt for an automated beta solver based on minimises residuals
% clear all;
% close all;

function [beta1, beta2, temp] = majorminor(indata, nfit, coords, res)
% coords = im_coords;
% nfit = n_approx;
data_size = size(indata);
npix = data_size(1)^2;
nside = data_size(1);
shapes = zeros(5,1);
midpt = round(nside/2);
radians = pi/180;
ncoeffs = 0.5*(nfit+1)*(nfit+2);
ang_res = res(1);
obj_size = res(2)*radians/60;

fcoeffs = zeros(ncoeffs, 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smooth/lower res data for faster running : aim is to capture most of the 
% in the first few basis functions

%nideal = 128;
%if nside > 200
%    data = resize(indata, nideal);
%else
    data = indata;
%end


k=0;
for i=1:nside
    for j=1:nside
        k=k+1;
        col_data(k,1)=data(i,j);
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coarse loop to isolate beta 1 and 2
imsig = 100;                                % init of MSE of current image
MSE = imsig;                                % MSE of current fit
prev_MSE = 2*imsig;                         % MSE of previous fit
prev_beta1 = 0.9*nfit^(-0.52)*obj_size;    % previous beta1
prev_beta2 = 0.9*nfit^(-0.52)*obj_size;    % previous beta2
max_beta = 0.9*nfit^(-0.52)*nside*ang_res;  % max allowed beta
beta2 = prev_beta2;                         % init of shapes
l=0;
step = 5*ang_res;                           % step size
change = 0;                                 % indicates convergence
c1=0;                                       % indicate change in beta1
c2=0;                                       % indicates change in beta2
m=1;

% tracks the optimisation routinue
temp(m,1)= prev_beta1;
temp(m,2)=prev_beta2;
temp(m,3)=MSE;

% coarse loop (step size of 5*ang_res) with loop limit of 10 iters
while (change < 2 && l<10)
    l=l+1;
    if c1 ~= 1;
	% initialise beta1 at guess + 5*step size
        beta1 = prev_beta1+5*step;

	% optimise beta1
        while (MSE < prev_MSE && m<1000) 
            m=m+1;
            prev_MSE = MSE;
            beta1 = beta1-step;
            if beta1 < 0
                beta1 = max_beta;
            end 
            MSE = shape_models(coords, nfit, beta1, beta2, col_data);
            temp(m,1)=beta1;
            temp(m,2)=beta2;
            temp(m,3)=MSE;
        end
        
	% is it a minimum? is it a global minimum?
        if round((beta1+step)*10e4) == round(prev_beta1*10e4)
            beta1=prev_beta1;
            c1=1;
        else
            prev_beta1 = beta1+step;
            beta1=prev_beta1;
            c2=0;
            c1=0;
        end
    end

    prev_MSE = 2*imsig;
    MSE = imsig;
    
    if c2 ~= 1
        beta2= prev_beta2+5*step; 

        % minimise beta2
        while (MSE < prev_MSE && m<1000)
            m=m+1;
            prev_MSE = MSE;
            beta2 = beta2-step; 
            if beta2 < 0
                beta2 = max_beta;
            end
            MSE = shape_models(coords, nfit, beta1, beta2, col_data);
            temp(m,1)=beta1;
            temp(m,2)=beta2;
            temp(m,3)=MSE;          
        end
    
	% is it a minimum? is it a global minimum?
        if round((beta2+step)*10e4) == round(prev_beta2*10e4)
            beta2=prev_beta2;
            c2 = 1;
        else
            prev_beta2 = beta2+step;
            beta2=prev_beta2;
            c2=0;
            c1=0;
        end 
    end

    change = c1+c2;
    prev_MSE = 2*imsig;
    MSE = imsig;    
end

m=m+1;
[sorts, ind] = sort(temp(:,3), 'ascend');
x = ind(1);
prev_beta1 = temp(x,1);
prev_beta2 = temp(x,2);

coarse=1
temp(m,1)=beta1;
temp(m,2)=beta2;
temp(m,3)=MSE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fine loop (step size = angular resolution)
imsig = 100;                                % init of MSE of current image
MSE = imsig;                                % MSE of current fit
prev_MSE = 2*imsig;   
step = ang_res;
force = 0; 				    % forces a beta value
c1=0;
c2=0;
change=0;

while (change < 2 && l<20)
    l=l+1;
    if c1 ~= 1;
        beta1 = prev_beta1+5*step;

        % minimise beta1
        while (MSE < prev_MSE && m<1000) 
            m=m+1;
            prev_MSE = MSE;
            beta1 = beta1-step;
            if beta1 < 0
                beta1 = beta1+step;
                force = 1;
            end
            MSE = shape_models(coords, nfit, beta1, beta2, col_data);
            temp(m,1)=beta1;
            temp(m,2)=beta2;
            temp(m,3)=MSE;
        end
       
        if round((beta1+step)*10e6) == round(prev_beta1*10e6)
            beta1=prev_beta1;
            c1=1;
        elseif force == 1
	        beta1=prev_beta1;
	        c1=1;
        else
            prev_beta1 = beta1+step;
            beta1=prev_beta1;
            c2=0;
            c1=0;
        end
    end

    prev_MSE = 2*imsig;
    MSE = imsig;
    
    force = 0;
    if c2 ~= 1
        beta2= prev_beta2+5*step; 

	% minimise beta2
        while (MSE < prev_MSE && m<1000)
            m=m+1;
            prev_MSE = MSE;
            beta2 = beta2-step; 
            if beta2 < 0
                beta2 = beta2+step;
		        force = 1;
            end         
            MSE = shape_models(coords, nfit, beta1, beta2, col_data);
            temp(m,1)=beta1;
            temp(m,2)=beta2;
            temp(m,3)=MSE;           
        end
    
        if round((beta2+step)*10e4) == round(prev_beta2*10e4)
            beta2=prev_beta2;
            c2 = 1;
        elseif force == 1
	        beta2=prev_beta2;
	        c2=1;
        else
            prev_beta2 = beta2+step;
            beta2=prev_beta2;
            c2=0;
            c1=0;
        end 
    end

    change = c1+c2;
    prev_MSE = 2*imsig;
    MSE = imsig;    
end

fine = 0
m=m+1;
temp(m,1)=beta1;
temp(m,2)=beta2;
temp(m,3)=MSE;

[sorts, ind] = sort(temp(:,3), 'ascend');
x = ind(1);
beta1 = temp(x,1);
beta2 = temp(x,2);


