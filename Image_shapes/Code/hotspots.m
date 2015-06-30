% Generate a matrix of parameters relating to a multi gaussian fit.

function [altered, gauss_fits, fail] = hotspots(altered, num_fits, ang_res)

data_size = size(altered);
xside = data_size(1);
yside = data_size(2);
micro = 0.2*xside;              % area used to fit gaussian
fail = 0;                       % was a fit achieved?

for g=1:num_fits
    new_max_pt = max(max(altered));
    for i=1:xside
        for j=1:yside
            if altered(i,j) == new_max_pt
                xmax = i;
                ymax = j;
            end
        end
    end
    
    k=0;
    nz = 0;
    small_gauss = zeros(micro, micro);          % area to fit
    small_coords = zeros(micro*micro,2);        % coords of area
    col_data = zeros(micro*micro,1);            % data as a column
    col_fit = zeros(micro*micro,1);             % fit as a column
    sq_fit = zeros(micro,micro);                % fit in 2D
    factor = 1;                                 % scale factor of fit
    
    % to create a smaller image (time saver)
    for i=1:micro
        for j=1:micro
            k=k+1;
            small_gauss(i,j) = altered(i+xmax-micro/2, j+ymax-micro/2);
            small_coords(k,1) = (i-micro/2)*ang_res;
            small_coords(k,2) = (j-micro/2)*ang_res;
            col_data(k,1) = small_gauss(i,j);
        end
    end
    
    small_gauss = small_gauss/new_max_pt;
    
    % Need to swap to a skewed gaussian fit
    for i=1:micro
        for j=1:micro 
             if small_gauss(i,j) <= 0.5
                small_gauss(i,j) = 0;
            else
                nz = nz+1;
            end
        end
    end   
    
    % Makes sure there is enough points before fitting
    if nz < 8
        fail = 1;
    else 
        [major, minor, PA, offsets] = cov_fit(small_coords, small_gauss);
    end
    
    paras = [big*sqrt(2*log(2)), small*sqrt(2*log(2)), -PA, offsets];
    
    col_fit = ell_gauss(paras, small_coords);
    sq_fit = convert(col_fit, micro, micro);

    % LS scale parameter on raw gaussian fit (not normalised)
    factor = inv(transpose(col_fit)*col_fit)*transpose(col_fit)*col_data;
    
    gauss_fits(i,1)=paras(1);
    gauss_fits(i,2)=paras(2);
    gauss_fits(i,3)=paras(3);
    gauss_fits(i,4)=(xmax+offsets(1))*ang_res;    
    gauss_fits(i,5)=(ymax+offsets(2))*ang_res;
    gauss_fits(i,6)=factor;
       
    for i=1:micro
        for j=1:micro 
            x = i+xmax-micro/2;
            y = j+ymax-micro/2;
            altered(x,y)= altered(x,y)-factor*sq_fit(i,j);
        end
    end
end