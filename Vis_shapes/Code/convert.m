% converts into a surface

function twoDsurf = convert( B, x_pxl, y_pxl)

twoDsurf = zeros(x_pxl, y_pxl);
k=1;
for x=1:x_pxl
    for y=1:y_pxl
        twoDsurf(x,y)=B(k);
        k=k+1;
    end
end

% figure(5)
% contour(real(twoDsurf));
% xlabel('x');
% ylabel('y');