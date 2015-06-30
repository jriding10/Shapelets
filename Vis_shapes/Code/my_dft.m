% takes visibility data and performs the 2D DFT -> works for the old data

function result = my_dft(im_coords, uv_coords, data)
% clear all;
% close all;
% 
% rawdata = load('PupA_vis_old.txt', '-ascii');
% shapes = load('PupA_dump_shapes.txt', '-ascii');
% im_coords = load('PupA_dump_imcoords.txt', '-ascii');
% uv_coords = load('PupA_dump_viscoords.txt', '-ascii');
% split_data = load('PupA_dump_fullmodel.txt', '-ascii');
% check = load('PupA_dump_image.txt','-ascii');
% 
% PA = shapes(3);
% uv_coords(:,1)=rawdata(:,1);
% uv_coords(:,2)=rawdata(:,2);
% PA=-pi/2;
% transform = [cos(PA), -sin(PA); sin(PA), cos(PA)];
% rot_coords = transform*transpose(uv_coords);
% uv_coords = transpose(rot_coords);
% uv_coords(:,1)=rawdata(:,1);
% uv_coords(:,2)=rawdata(:,2);
% radians = pi/180;
% im_res = 0.001;
% pxl = 128;
% 
% k=0;
% for i=1:pxl
%     for j=1:pxl
%         k=k+1;
%         im_coords(k,1)=((i-1)-pxl/2)*im_res;
%         im_coords(k,2)=((j-1)-pxl/2)*im_res;
%     end
% end

% im_res = 0.001;

uvlen = size(uv_coords);
imlen = size(im_coords);
pxl = sqrt(imlen(1));

% PA=shapes(3);
% transform = [cos(PA), -sin(PA); sin(PA), cos(PA)];
% rot_coords = transform*transpose(uv_coords);
% uv_coords = transpose(rot_coords);

% for i=1:uvlen(1)
%     data(i)=rawdata(i,4)+1i*rawdata(i,5);
% end

result = zeros(pxl);
mat_result = zeros(imlen(1));

for i=1:uvlen(1)
    for j=1:imlen(1)
        p_term = uv_coords(i,1)*im_coords(j,1)+uv_coords(i,2)*im_coords(j,2);
        mat_result(j)=mat_result(j)+1/uvlen(1)*data(i)*exp(1i*2*pi*p_term);
    end
end

k=0;
for i=1:pxl
    for j=1:pxl
        k=k+1;
        result(i,j)=mat_result(k);
    end
end

% figure(1)
% contour(result);
% figure(2)
% contour(check);
% res = check-result;
% figure(3)
% surf(abs(res));
        