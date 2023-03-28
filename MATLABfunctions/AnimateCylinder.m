function AnimateCylinder(VORT,varargin)
%takes stacked datamatricies [X,Y,t] and animates throuth time t. 
%varargin {1} defines the timestep used when saving animation. 
%any varargin {2} initiates a gif saving.  

f1 = figure;

vortmin = -5;  % only plot what is in -5 to 5 range
vortmax = 5;
VORT(VORT>vortmax) = vortmax;  % cutoff at vortmax
VORT(VORT<vortmin) = vortmin;  % cutoff at vortmin

plot_cylinder = imagesc(VORT(:,:,1)); % plot vorticity field
load CCcool.mat 
colormap(CC);  % use custom colormap


% clean up axes
set(gca,'XTick',[1 50 100 150 200 250 300 350 400 449],'XTickLabel',{'-1','0','1','2','3','4','5','6','7','8'})
set(gca,'YTick',[1 50 100 150 199],'YTickLabel',{'2','1','0','-1','-2'});

axis equal
hold on

% add contour lines (positive = solid, negative = dotted)
[~,con1] = contour(VORT(:,:,1),[-5.5:.5:-.5 -.25 -.125],':k','LineWidth',1.2);
[~,con2] = contour(VORT(:,:,1),[.125 .25 .5:.5:5.5],'-k','LineWidth',1.2);

theta = (1:100)/100'*2*pi;
x = 49+25*sin(theta);
y = 99+25*cos(theta);
fill(x,y,[.3 .3 .3])  % place cylinder
plot(x,y,'k','LineWidth',1.2) % cylinder boundary

set(gcf,'PaperPositionMode','auto') 
xlabel('x'); ylabel('y'); % ax labels


for i = 1: size(VORT,3)
    plot_cylinder.CData = VORT(:,:,i);
    con1.ZData = VORT(:,:,i);
    con2.ZData = VORT(:,:,i);

    pause(0.005) 

     % make .gif animation
    if nargin == 3
        frame = getframe(gcf); % gcf = get current figure
        image = frame2im(frame);
        filename = 'Cylinder_flow_anim.gif';
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
