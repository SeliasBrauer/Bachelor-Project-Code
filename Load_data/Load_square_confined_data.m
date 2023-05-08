% import and animate square flow: 
clc; clear variables; close all; 
dx = 1/30 ;
string1 = 'C:\Users\selia\Desktop\Simulation data';
time = 5001 : 1 : 7500;

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
%%
%AnimateSquare(vg,'Square_confined_anim.gif')
VORTALL_CONFINED_SINDy_Large = reshape(vg,[],length(time)) ;

save("VORTALL_CONFINED_SINDy_Large","VORTALL_CONFINED_SINDy_Large");
