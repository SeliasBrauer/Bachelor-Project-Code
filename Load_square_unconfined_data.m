% import and animate square flow: 
clc; clear variables; close all; 
dx = 1/22 ;
string1 = 'C:\Users\selia\Desktop\vorticity unconfined-20230303T215318Z-001\vorticity unconfined\vorticity_table_';
time = 5010 : 10 : 7500;

for i = 1:length(time)
    
    string2 = sprintf('%i.csv',time(i));
    path = append(string1,string2);
    Data = readmatrix(path);
    Vorticity = Data(:,1); 
    x = Data(:,2);
    y = Data(:,3);

    xg = min(x) : dx : max(x); 
    yg = flip( (min(y) : dx : max(y))', 1); 

    [Xg,Yg,vg(:,:,i)] = griddata(x,y,Vorticity,xg,yg);
end


%set vorticity within cylinder boundary to zero
for i = 1:length(xg)
    for j = 1:length(yg)
        if (xg(i) > 0 && xg(i) < 1) && (yg(j) < 0.5 && yg(j) > -0.5)
            vg(j,i,:) = 0;
        end 
    end
end 
save("VORSTACK","vg");