function AnimateSquare(VORT,filename,varargin)
%takes stacked datamatricies [X,Y,t] and animates throuth time t. 
%varargin {1} defines the timestep used when saving animation. 
%any varargin {2} initiates a gif saving. 

%domain properties
xstart = -1; 
ystart = 2; 
dx = 1/22;

f1 = figure(); 
vortmin = -5;  % only plot what is in -5 to 5 range
vortmax = 5;
VORT(VORT>vortmax) = vortmax;  % cutoff at vortmax
VORT(VORT<vortmin) = vortmin;  % cutoff at vortmin

plot_square = imagesc(VORT(:,:,1));
load CCcool.mat 
colormap(CC);  % use custom colormap
hold on; axis equal; 
colorbar(); 

% add contour lines (positive = solid, negative = dotted)
[~,con1] = contour(VORT(:,:,1),[-5.5:.5:-.5 -.25 -.125],':k','LineWidth',1.2);
[~,con2] = contour(VORT(:,:,1),[.125 .25 .5:.5:5.5],'-k','LineWidth',1.2);


% clean up axes
Xtick = [1 : 1/dx : size(VORT(:,:,1),2) ]; 
Ytick = [1 : 1/dx : size(VORT(:,:,1),1) ]; 

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


%plot square figure: 
x1 = -xstart*1/dx + 1 ; x2 = x1 + 1/dx;
x = [x1, x2, x2, x1, x1];

y1 = (ystart - 0.5)*1/dx + 1;  y2 = y1 + 1/dx;
y = [y1, y1, y2, y2, y1 ];

fill(x,y,[.3 .3 .3])  % place square
plot(x,y,'k','LineWidth',1.2) % square boundary

xlabel('x'); ylabel('y'); % ax labels


for i = 1: size(VORT,3)
    plot_square.CData = VORT(:,:,i);
    con1.ZData = VORT(:,:,i);
    con2.ZData = VORT(:,:,i);

    pause(0.005) 

     % make .gif animation
    if nargin == 4
        frame = getframe(gcf); % gcf = get current figure
        image = frame2im(frame);
        %filename = 'Square_flow_anim.gif';
        dt = varargin{1};
        [A,map] = rgb2ind(image,256);
        if i == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',dt,'BackgroundColor',0);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',dt);
        end
    end

end

end