% Noise estimator: takes the area including the object, then deletes the
% objects and places non-object data into a vector. Needs to be told start
% and stop values - ie not smart

% % static values for Fornax A
% nsize = 300;                % size of noise image
% r_center = 500;             % row center of data
% c_center = 500;             % column center of data
% r_imstart = 120;            % row start of image
% r_imstop = 220;             % row end of image
% c_imstart = 100;            % column start of image
% c_imstop = 230;             % column end of image
% nsum = 0;
% norm = 0;
% nim = zeros(nsize);

% static values for Puppis A
nsize = 100;                % size of noise image
r_center = 150;             % row center of data
c_center = 150;             % column center of data
r_imstart = 20;             % row start of image
r_imstop = 80;             % row end of image
c_imstart = 20;            % column start of image
c_imstop = 75             % column end of image
nsum = 0;
norm = 0;
nim = zeros(nsize);
noise = zeros(nsize);

% % static values for Pup A dump
% nsize = 64;                % size of noise image
% r_center = 64;             % row center of data
% c_center = 64;             % column center of data
% r_imstart = 24;            % row start of image
% r_imstop = 40;             % row end of image
% c_imstart = 32;            % column start of image
% c_imstop = 44;             % column end of image
% nsum = 0;
% norm = 0;
% nim = zeros(nsize);

rstart = r_center - nsize/2;
cstart = c_center - nsize/2;

for i=1:nsize
    for j=1:nsize
        noise(i,j)=mid_im(i+rstart, j+cstart);
    end
end

k=0;
for i=1:nsize
    for j=1:nsize
        if j<c_imstart
            k=k+1;
            nsum=nsum+noise(i,j)^2;
            norm=norm+noise(i,j);
            nim(i,j)=noise(i,j);
            nvec(k)=noise(i,j);
        end
        if (i>r_imstop && j>=c_imstart && j<=c_imstop)
            k=k+1;
            nsum = nsum+noise(i,j)^2;
            norm=norm+noise(i,j);
            nim(i,j)=noise(i,j);
            nvec(k)=noise(i,j);
        end
        if (i<r_imstart && j>=c_imstart && j<=c_imstop)
            k=k+1;
            nsum = nsum+noise(i,j)^2;
            norm=norm+noise(i,j);
            nim(i,j)=noise(i,j);
            nvec(k)=noise(i,j);
        end
        if j>c_imstop
            k=k+1;
            nsum = nsum+noise(i,j)^2;
            norm=norm+noise(i,j);
            nim(i,j)=noise(i,j);
            nvec(k)=noise(i,j);
        end
    end
end

% check normalisation
% nmse = var(nvec)/max(max(mid_im));
nmse = var(nvec)
nmean = mean(nvec)