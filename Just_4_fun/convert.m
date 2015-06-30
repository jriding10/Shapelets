% converts coords and data into xaxis, yaxis and immage

function immage = convert(data, x, y)

immage = zeros(x,y);
xaxis = zeros(x,1);
yaxis = zeros(y,1);

k=0;
for i=1:x
    for j=1:y
        k=k+1;
        immage(i,j)=data(k);
    end
end
