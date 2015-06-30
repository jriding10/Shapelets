% Takes the miriad + excel output and creates a useable matlab file File
% format: vis #, time (fraction of day in seconds), ant, ant, u, v, amp,
% phase (8) but read in as a single column. There are 32 baselines + self
% cal => 528 measurements every 2 seconds

% To do: create rows and columns, delete flagged data, integrate up to 8
% seconds
clear all;
close all;

rawdata = load('../Common/PupA.txt', '-ascii');

s = size(rawdata);
rows = s(2)/8;
inttime = 2;
rts_time = 8;
diff_int = rts_time/inttime;
nxt_vis = 528;                          % number of baselines

better_data = zeros(rows, 8);           % raw data as it is suppose to be 

% columnate the data sensibly
k=0;
flagged=0;
for i=1:rows
    for j=1:8
        k=k+1;
        better_data(i,j)=rawdata(k);
    end
    if better_data(i,7)==0
         flagged=flagged+1;
    end    
end

% remove flagged data, express time as time from the start of obs, express
% reading as complex
useful_data = zeros(rows-flagged, 8);     

rows2=0;
for i=1:rows
    if (better_data(i,7) ~= 0) || (better_data(i,8) ~= 0)
        rows2=rows2+1;        
        useful_data(rows2,1) = better_data(i,1);
        useful_data(rows2,2) = round((better_data(i,2)-better_data(1,2))*24*60*60);
        useful_data(rows2,3) = abs(better_data(i,3));
        useful_data(rows2,4) = better_data(i,4);
        useful_data(rows2,5) = better_data(i,5);
        useful_data(rows2,6) = better_data(i,6);
        useful_data(rows2,7) = better_data(i,7)*exp(1i*better_data(i,8));
    end
end

start_time=useful_data(1,2)+better_data(1,2);     % percent seconds of day

% integrate 8 second intervals: atm remove extra at the end
rows3 = floor(rows2/diff_int);
final_data = zeros(rows3, 6);
row_set = zeros(rows3,5);           % testing

% ant 1 is flagged as are self cal, so max number of baselines = 465. To
% find approx number based on first measurement set;

unique=1;
while useful_data(unique,2) == useful_data(1,2)
    unique=unique+1;
end
nxt_row = unique;

row_num = 0;
% only have to check vis # and assumes that baselines do not move over 8
% second intervals
for i=1:rows2
    if useful_data(i,8)==0
        row_num=row_num+1;
        vis_num = useful_data(i,1);
        sum = useful_data(i,7);
        score = 1;
        row_set(row_num, 2) = row_num;
        for j=1:3
            start = i+j*nxt_row-32;
            stop = i+j*nxt_row+32;
            thing = 2;
            for k = start:stop
                if k < rows2
                    if useful_data(k,1) == vis_num+j*528
                        thing = thing+1;
                        sum=sum+useful_data(k,7);
                        useful_data(k,8)=1;
                        score = score+1;
                        row_set(row_num,thing)=k;
                    end
                end
            end
        end
        row_set(row_num,1)=score;
        temp=sum/score;
        final_data(row_num,1)=useful_data(i,5)*1000;
        final_data(row_num,2)=useful_data(i,6)*1000;
        final_data(row_num,4)=real(temp);
        final_data(row_num,5)=imag(temp);
        final_data(row_num,3)=sqrt(final_data(row_num,1)^2+final_data(row_num,2)^2);
        final_data(row_num,6)=score/4;
    end
end
  
% outputs the first 8 seconds
snapshot = zeros(unique, 6);
for i=1:unique
    snapshot(i,:)=final_data(i,:);
end

save('PupA_vis_full.txt', 'final_data', '-ascii', '-double');
save('PupA_vis_snapshot.txt', 'snapshot', '-ascii', '-double');
                    
                        
                
        
