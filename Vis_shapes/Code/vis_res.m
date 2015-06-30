% Old file that computed the IDFT for dump data

% function [thing1, thing2, thing3] = vis_res(coords, imCoords, shapes, fout, nmax, weight, data) 
clear all;
close all;

coords = load('../Text/PupA_dump_coords.txt', '-ascii');
shapes = load('../Text/PupA_shapes_flipped.txt', '-ascii');
fcoeffs = load('../Text/PupA_coeffs_flipped.txt', '-ascii');
weight = load('../Text/PupA_dump_weights.txt', '-ascii');
uvdata = load('../Text/PupA_dump_recodata.txt', '-ascii');
check = load('../Text/PupA_dump_image.txt', '-ascii');

n1max = max(fcoeffs(:,1));
n2max = max(fcoeffs(:,2));
nmax = max(n1max, n2max);
relcoeffs = 200;

uvlen = size(coords);
pxl = 128;
tot_pxl = pxl*pxl;
coeffs = size(fcoeffs);
radians = pi/180;
ang_res = 0.078125*radians;
FoV =  10*radians;
% shapes = ang_res*shapes;
shapes(4) = -1*ang_res;
shapes(5) = -1*ang_res;

k=0;
for i=1:pxl
    for j=1:pxl
        k=k+1;
        imCoords(k,1)=-FoV/2+i*ang_res;
        imCoords(k,2)=-FoV/2+j*ang_res;
    end
end

uvmodel = zeros(uvlen(1), 1);
I = zeros(tot_pxl,1);
I_m = zeros(tot_pxl,1);
I_r = zeros(tot_pxl,1);

type = 1;
M = basis(coords, shapes, type, nmax);

[sorts, ind] = sort(abs(fcoeffs(:,3)), 'descend');

for i=1:relcoeffs
    num = ind(i);
    fout(i,1) = fcoeffs(num,1);
    fout(i,2) = fcoeffs(num,2);
    fout(i,3) = fcoeffs(num,3);
end 

for i=1:relcoeffs(1)
    n1 = fout(i,1);
    n2 = fout(i,2);
    f_hat = fout(i,3);
    Mpiece = M(:,n1+1,n2+1);
    uvmodel = uvmodel+f_hat*(Mpiece);
end

weight = ones(uvlen(1),1);

tot_w = sum(weight);
for i=1:uvlen(1)
    data(i) = uvdata(i,1)+1i*uvdata(i,2);
end

% fudge = (2.5623e+05);
fudge = 1;
fudge = 2.6760e+06 - 3.2100e+05i;

Imat = 1/tot_w*data.'.*(weight);
Imat_m = fudge/tot_w*(0.5*uvmodel).*(weight);
Imat_r = Imat - Imat_m;

for i=1:tot_pxl
   for j=1:uvlen(1)
       I(i)   = I(i)  +Imat(j)  *exp(2*pi*1i*(coords(j,1)*imCoords(i,1)+coords(j,2)*imCoords(i,2)));
       I_m(i) = I_m(i)+Imat_m(j)*exp(2*pi*1i*(coords(j,1)*imCoords(i,1)+coords(j,2)*imCoords(i,2)));
       I_r(i) = I_r(i)+Imat_r(j)*exp(2*pi*1i*(coords(j,1)*imCoords(i,1)+coords(j,2)*imCoords(i,2)));

   end
end

% I_m = 355*I_m;

data_image =zeros(pxl);
model_image = zeros(pxl);
resid_image = zeros(pxl);
xyaxis = zeros(pxl,1);

k=0;
for i=1:pxl  
    for j=1:pxl
        k=k+1;
        data_image(i,j)=real(I(k));
        model_image(i,j)=real(I_m(k));
        resid_image(i,j)=real(I_r(k));
    end
end

for i=1:pxl
    xyaxis(i) = imCoords(i,2)/radians;
end

% figure(1)
% contour(check);
% title('Answer');
% 
% figure(2);
% contour(data_image)
% title('UV data to Image');
% hold on;
% plot(64, 64, '*');
% hold off;
% 
% figure(3);
% contour(model_image)
% hold on;
% title('UV Model to Image');
% plot(64, 64, '*');
% hold off;
% 
% figure(4);
% contour(resid_image)
% title('UV Resids to Image');

for i=56:72
    for j=50:66
        new_data(i-55, j-49) = data_image(i,j);
        new_model(i-55, j-49) = model_image(i,j);
    end
end

[NMSE, PSNR, SSIM] = stat_calc(new_data, new_model)

% new_size = 64;
% mid = new_size/2;
% for i=1:new_size
%     naxis(i)=xyaxis(i+mid-1);
%     for j=1:new_size
%         ndata(i,j)=data_image(i+mid-1, j+mid-1);
%         nmodel(i,j)=model_image(i+mid-1, j+mid-1);
%         nresids(i,j)=resid_image(i+mid-1, j+mid-1);
%     end
% end

cstart = 47;
cstop = 72;
rstart = 50;
rstop = 75;
diff = cstop-cstart;
czero = 64-cstart;
rzero = 64-rstart;

for i=rstart:rstop
    xaxis(i-rstart+1)=((i-rstart+1)-rzero)*0.078125;
    yaxis(i-cstart-2)=-((i-cstart+1)-czero)*0.078125;
    for j=cstart:cstop
        ndata(i-rstart+1, j-cstart+1) = data_image(i,j);
        nmodel(i-rstart+1, j-cstart+1)= model_image(i,j);
        nresids(i-rstart+1, j-cstart+1)= resid_image(i,j);
    end
end
    
figure(1);
surf(yaxis, xaxis, ndata, 'LineStyle','none');
view(0,90);
axis image;
xlabel('Degrees','FontSize',30);
ylabel('Degrees','FontSize',30);

figure(2);
surf(yaxis, xaxis, nmodel, 'LineStyle','none');
view(0,90);
axis image;
xlabel('Degrees','FontSize',30);
ylabel('Degrees','FontSize',30);

figure(3);
surf(yaxis, xaxis, nresids, 'LineStyle','none');
view(0,90);
axis image;
xlabel('Degrees','FontSize',30);
ylabel('Degrees','FontSize',30);





% figure(1);
% contour(ndata);
% hold on;
% plot(13,13,'*');
% hold off;
% 
% figure(2);
% contour(nmodel);
% hold on;
% plot(13,13,'*');
% hold off;
% 
% figure(3);
% contour(nresids);
% hold on;
% plot(13,13,'*');
% hold off;