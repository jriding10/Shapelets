% Deal with the RTS dump data to create a compadible txt file

% col: u, v, weight, real, imag
data = load('../Text/PupA_dump.txt', '-ascii');
s = size(data);
l = 1;
minbl = 30;

% count the number of unique measurements
unique = 1;
while l ~= s(1)
    if data(l+1,1) == data(l,1) && data(l+1,2) == data(l,2)
        l=l+1;
    else
        l=l+1;
        unique = unique+1;
    end
end

% calcs bl and integrates multiple measurements
avg_vis = zeros(unique,6);	 
l = 1;
k = 1;
j = 1;

while l ~= s(1)
    if data(l+1,1) == data(l,1) && data(l+1,2) == data(l,2) 
        avg_vis(j,3) = avg_vis(j,3)+data(l,3);        
        avg_vis(j,4) = avg_vis(j,4)+data(l,4);
        avg_vis(j,5) = avg_vis(j,5)+data(l,5);
        l=l+1;
        k=k+1;
        
    else 
        avg_vis(j,1) = data(l,1);
        avg_vis(j,2) = data(l,2);
	    avg_vis(j,6) = sqrt(power(avg_vis(j,1),2)+power(avg_vis(j,2),2));
        avg_vis(j,3) = avg_vis(j,3)+data(l,3);
        avg_vis(j,4) = (avg_vis(j,4)+data(l,4))/k;
        avg_vis(j,5) = (avg_vis(j,5)+data(l,5))/k;
        
        l=l+1;
        k=1;
        j=j+1;
    end

end

% fix to get the last row right
avg_vis(j,1) = data(l,1);
avg_vis(j,2) = data(l,2);
avg_vis(j,6) = sqrt(power(avg_vis(j,1),2)+power(avg_vis(j,2),2));
avg_vis(j,3) = avg_vis(j,3)+data(l,3);
avg_vis(j,4) = (avg_vis(j,4)+data(l,4))/k;
avg_vis(j,5) = (avg_vis(j,5)+data(l,5))/k;

norm = max(avg_vis(:,3));

% turn data into a more convenient/readable form	   
final_data = zeros(unique,6);

% rearrange data for use with Vis_shapelets
final_data(:,1) = avg_vis(:,1);
final_data(:,2) = avg_vis(:,2);
final_data(:,3) = avg_vis(:,6);
final_data(:,4) = avg_vis(:,4);
final_data(:,5) = avg_vis(:,5);
final_data(:,6) = avg_vis(:,3)/norm;

save('PupA_vis_old.txt', 'final_data', '-ascii');
