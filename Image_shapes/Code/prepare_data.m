% Prepare data: First pass at data - determines side lengths, zero pads,
% zooms and creates first pxl and ra/dec coordinate axes

function [coords, outdata, col_data, sky_coords] = prepare_data(pxlobj, obj, ares, indata)

%indata = data;
%ares = ang_res;
%pxlobj = pix_obj;
%obj = [coord_obj, obj_size_max];

cutoff = 0.05;
radians = pi/180;
fudge = 15;

% Sizings
data_size = size(indata);
max_flux = max(max(indata));
obj_size = round(((obj(3)/60*radians)/ares)*fudge-1); 
npix_side = min(data_size(1), data_size(2));
ra_pxl = pxlobj(1);
dec_pxl = pxlobj(2);

% Should we zoom? - Image is much greater than source (and position)
if obj_size > npix_side
    nside = npix_side;
else
    nside = obj_size;
end

% makes sure the resized dataset falls within original data
midpix = floor((nside+1)/2);
midpix_ra = 0;
midpix_dec = 0;

if (ra_pxl+midpix)>npix_side
    midpix_ra = npix_side-ra_pxl;
end
if (ra_pxl-midpix)<1
    midpix_ra = ra_pxl;
end
if midpix_ra == 0
    midpix_ra = midpix;
end
nside_ra = 2*midpix_ra;

if (dec_pxl+midpix)>npix_side
    midpix_dec = npix_side-dec_pxl;
end
if (dec_pxl-midpix)<1
    midpix_dec = dec_pxl;
end
if midpix_dec == 0
    midpix_dec = dec_pxl;
end
nside_dec = 2*midpix_dec;

nside = min(nside_ra, nside_dec);
midpix = nside/2; 
% columns are ra so col_axis pertains to matrix position, not label
% position

for i=1:nside
    rowaxis(i,1) = (i-midpix);
    colaxis(i,1) = (i-midpix);
end

k=0;
for i=1:nside
    for j=1:nside     
        k=k+1;
        coords(k,1)=rowaxis(i)*ares;
        coords(k,2)=colaxis(j)*ares;
        outdata(i,j) = indata(i+pxlobj(2)-midpix, j+pxlobj(1)-midpix);
	%if outdata(i,j) < cutoff*max_flux
        %    outdata(i,j)=0;
        %end
        col_data(k,1) = outdata(i,j);    
    end
end

%dec = row_axis/radians + obj(2)*ones(nside,1);
%ra = -col_axis/radians - obj(1)*ones(nside,1);

dec = rowaxis;
ra = colaxis;
sky_coords = [ra, dec];
