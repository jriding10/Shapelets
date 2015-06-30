% Retrieves information about the fits file

function [pix_obj, coord_obj, ang_res]=fits_info(info)
% lots of info in a fits file, starting at line 10 should be fine
% (CRPIXx, CDELTx, CRVALx) = (CR=control)(pixel number, resolution, pixel
% value)(1=RA, 2=Dec).
radians = pi/180;
sinfo = size(info.PrimaryData.Keywords);
max_search = sinfo(1);
all = 0;

x=1;
while all ~= 6
    temp = info.PrimaryData.Keywords(x);
    CP1 = strcmp(temp, 'CRPIX1');
    CV1 = strcmp(temp, 'CRVAL1');
    CD1 = strcmp(temp, 'CDELT1') | strcmp(temp, 'CD1_1');
    if CP1 ==1
        CRPIX1 = info.PrimaryData.Keywords{x,2};
	all = all+1;
    end
    if CD1 ==1
        CDELT1 = info.PrimaryData.Keywords{x,2};
	all = all+1;
    end
    if CV1 ==1
        CRVAL1 = info.PrimaryData.Keywords{x,2};
	all = all+1;
    end
    CP2 = strcmp(temp, 'CRPIX2');
    CV2 = strcmp(temp, 'CRVAL2');
    CD2 = strcmp(temp, 'CDELT2') | strcmp(temp, 'CD2_2');
    if CP2 ==1
        CRPIX2 = info.PrimaryData.Keywords{x,2};
	all = all+1;
    end
    if CD2 ==1
        CDELT2 = info.PrimaryData.Keywords{x,2};
	all = all+1;
    end
    if CV2 ==1
        CRVAL2 = info.PrimaryData.Keywords{x,2};
	all = all+1;
    end
    x = x+1;
    if x > max_search && all < 6
	error('Failure to obtain fits values');
    end
end 

yres = abs(CDELT1);
xres = abs(CDELT2);

% Locate object in terms of pixel number from ra/dec info
ang_res = yres*radians;             % assumes yres == xres
pix_obj = [CRPIX1,CRPIX2];          % [ra, dec] as pixel number 
coord_obj = [CRVAL1, CRVAL2];       % [ra, dec] as degrees
