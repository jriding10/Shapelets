function [beta1, beta2] = circ_beta(col_data, nfit, coords, res)
data_size = size(col_data);
npix = data_size(1);
nside = sqrt(npix);
shapes = zeros(5,1);
midpt = round(nside/2);
radians = pi/180;
ncoeffs = 0.5*(nfit+1)*(nfit+2);
ang_res = res(1);
obj_size = res(2)*radians/60;

fcoeffs = zeros(ncoeffs, 3);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSE = 100;                                % MSE of current fit
prev_MSE = 2*MSE;                         % MSE of previous fit
beta = 0.9*nfit^(-0.52)*obj_size;        % previous beta
max_beta = 0.9*nfit^(-0.52)*nside*ang_res;
l=1;
step = 5*ang_res;                           % step size
beta = beta+5*step;

while (MSE < prev_MSE && l<1000)
    prev_MSE = MSE;
    beta = beta-step;
    if (beta < 0)
        beta = max_beta;
    end    
    MSE = shape_models(coords, nfit, beta, beta, col_data);
    result(l,1) = beta;
    result(l,2) = MSE;
    result(l,3) = prev_MSE;
    l=l+1;
end

step = ang_res;                           % step size
beta = beta+5*step;

while (MSE < prev_MSE && l<1000)
    prev_MSE = MSE;
    beta = beta-step;
    if (beta < 0)
        beta = max_beta;
    end    
    MSE = shape_models(coords, nfit, beta, beta, col_data);
    result(l,1) = beta;
    result(l,2) = MSE;
    result(l,3) = prev_MSE;
    l=l+1;
end

result(l,1) = beta;
result(l,2) = MSE;
result(l,3) = prev_MSE;
beta1=result(l-2,1);
beta2=beta1;
