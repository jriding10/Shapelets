% takes visibility data and performs the 2D DFT -> works for the old data

function result = my_dft(imCoords, coords, data)
uvlen = size(coords);
imlen = size(imCoords);
pxl = sqrt(imlen(1));
I = zeros(imlen(1),1);

% PA=shapes(3);
% transform = [cos(PA), -sin(PA); sin(PA), cos(PA)];
% rot_coords = transform*transpose(uv_coords);
% uv_coords = transpose(rot_coords);

% for i=1:uvlen(1)
%     data(i)=rawdata(i,4)+1i*rawdata(i,5);
% end

result = zeros(pxl);
mat_result = zeros(imlen(1));

% for i=1:uvlen(1)
%     for j=1:imlen(1)
%         p_term = uv_coords(i,1)*im_coords(j,1)+uv_coords(i,2)*im_coords(j,2);
%         mat_result(j)=mat_result(j)+1/uvlen(1)*data(i)*exp(1i*2*pi*p_term);
%     end
% end
tmp_wgts = ones(uvlen(1),1);

for j=1:imlen(1)
   for k=1:uvlen(1)
       I(j)=I(j)+1/sum(tmp_wgts)*data(k)*tmp_wgts(k)*exp(2*pi*1i*(coords(k,1)*imCoords(j,1)+coords(k,2)*imCoords(j,2)));
   end
end


k=0;
for i=1:pxl
    for j=1:pxl
        k=k+1;
        result(i,j)=I(k);
    end
end

% figure(1)
% contour(result);
% figure(2)
% contour(check);
% res = check-result;
% figure(3)
% surf(abs(res));
        