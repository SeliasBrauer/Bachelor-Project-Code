function f1 = plotSquare(VORT,dx)
%domain properties
xstart = -1; 
ystart = 2; 

f1 = figure(); 
vortmin = -5;  % only plot what is in -5 to 5 range
vortmax = 5;
VORT(VORT>vortmax) = vortmax;  % cutoff at vortmax
VORT(VORT<vortmin) = vortmin;  % cutoff at vortmin

imagesc(VORT)
load CCcool.mat 
colormap(CC);  % use custom colormap
hold on; axis equal; 
colorbar(); 

% add contour lines (positive = solid, negative = dotted)
contour(VORT,[-5.5:.5:-.5 -.25 -.125],':k','LineWidth',1.2)
contour(VORT,[.125 .25 .5:.5:5.5],'-k','LineWidth',1.2)


% clean up axes
Xtick = [1 : 1/dx : size(VORT,2) ]; 
Ytick = [1 : 1/dx : size(VORT,1) ]; 

for i = 1:length(Xtick)
    x = xstart + i - 1;
    Xtick_label{i} = sprintf('%i', x);
end

for i = 1:length(Ytick)
    y = ystart + 1 - i;
    Ytick_label{i} = sprintf('%i', y);
end


set(gca,'XTick',Xtick,'XTickLabel',Xtick_label);
set(gca,'YTick',Ytick,'YTickLabel',Ytick_label);



x1 = -xstart*1/dx + 1 ; x2 = x1 + 1/dx;
x = [x1, x2, x2, x1, x1];

y1 = (ystart - 0.5)*1/dx + 1;  y2 = y1 + 1/dx;
y = [y1, y1, y2, y2, y1 ];

fill(x,y,[.3 .3 .3])  % place square
plot(x,y,'k','LineWidth',1.2) % square boundary

xlabel('x'); ylabel('y'); % ax labels

end
