% converts my line data into an image

function out = views(data, xbit)

len = size(data);
pxl = sqrt(len(1));

square = convert(data, pxl, pxl);
for i=1:xbit
    for j=1:xbit
        out(i,j)=square((pxl-xbit)/2+i, (pxl-xbit)/2+j);
    end
end

figure(10)
surf(out);
view(0,90);