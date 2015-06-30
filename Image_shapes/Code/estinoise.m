% estimate the rms noise near the source

function [nmse, check] = estinoise(data)

data_size = size(data);
nside = data_size(1);

edge_margin = 0.1;
delta = 0.1*nside;

ntot = 0;

k=0;
for i=1:nside
    for j=1:nside
       if i < delta
           k=k+1;
           noise(k)=data(i,j);
           cutout(i,j)=data(i,j);
           ntot = ntot+noise(k);
        end
        if i > (nside-delta)
            k=k+1;
            noise(k)=data(i,j);
            cutout(i,j)=data(i,j);
            ntot = ntot+noise(k);
        end
        if j < delta && i > delta && i < (nside - delta)
            k=k+1;
            noise(k)=data(i,j);
            cutout(i,j)=data(i,j);
            ntot = ntot+noise(k);
        end
        if j > (nside-delta) && i > delta && i < (nside - delta)
            k=k+1;
            noise(k)=data(i,j);
            cutout(i,j)=data(i,j);
            ntot = ntot+noise(k);
        end
    end
end

nmse = var(noise);
check = (ntot-mean(noise)^2)/k;
