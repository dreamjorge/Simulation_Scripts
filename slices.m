% slices aproach
znorm = z;
xnorm = x;
ynorm = y;

xmin   = min(znorm(:));
xmax   = max(znorm(:)); 
ymin   = min(xnorm(:));  
ymax   = max(xnorm(:)); 
zmin   = min(xnorm(:));  
zmax   = max(xnorm(:));

bb         = Wo;
 iindex     = find(abs(bb)<=.0005);
 bb(iindex) = NaN;

% hslice = surf(linspace(xmin,xmax,100),linspace(ymin,ymax,100),zeros(100));
% rotate(hslice,[-1,0,0],-45)
% 
% xd1     = get(hslice,'XData');
% yd1     = get(hslice,'YData');
% zd1     = get(hslice,'ZData');
% 
% delete(hslice)
% 
% hslice = surf(linspace(xmin,xmax,100),linspace(ymin,ymax,100),zeros(100));
% rotate(hslice,[-1,0,0],45)


% xd2     = get(hslice,'XData');
% yd2     = get(hslice,'YData');
% zd2     = get(hslice,'ZData');
% 
% delete(hslice)
close(figure(101))
figure (101)
%plotPropagatedRays(rayH11,rayH12,rayH21,rayH22)
%slice to 45 degrees
% h45 = slice(z,x,y,abs(Wo),xd1,yd1,zd1); 
% h45.FaceColor = 'interp'; 
% h45.EdgeColor = 'none'; 
% h45.DiffuseStrength = 0.15;

plotPropagatedRays(rayH11,rayH12,rayH21,rayH22,1,1)
hold on 
%slice to -45 degrees
% h_45 = slice(z,x,y,abs(Wo),xd2,yd2,zd2); 
% h_45.FaceColor = 'interp'; 
% h_45.EdgeColor = 'none'; 
% h_45.DiffuseStrength = 0.15;

%Border of slices
% hxmin = slice(znorm,xnorm,ynorm,abs(bb),xmin,[],[]); 
% hxmin.FaceColor = 'interp'; 
% hxmin.EdgeColor = 'none'; 
% hxmax = slice(znorm,xnorm,ynorm,abs(bb),xmax,[],[]); 
% hxmax.FaceColor = 'interp'; 
% hxmax.EdgeColor = 'none';
% hy = slice(znorm,xnorm,ynorm,abs(bb),[],ymax,[]); 
% hy.FaceColor = 'interp'; 
% hy.EdgeColor = 'none';  
% 
% hz = slice(znorm,xnorm,ynorm,abs(bb),[],[],zmin); 
% hz.FaceColor = 'interp'; 
% hz.EdgeColor = 'none';

%% obstruction 
% xobs = [-0.5*lx,0,0.5*lx]; 
% yobs = [-0.5*ly,0.5,0.5*ly];
% 
% hxobs = slice(z,x,y,abs(Wo),[],xobs,[]); 
% hxobs(1).FaceColor = 'interp'; 
% hxobs(2).FaceColor = 'interp'; 
% hxobs(3).FaceColor = 'interp'; 
% hxobs(1).EdgeColor = 'none'; 
% hxobs(2).EdgeColor = 'none';
% hxobs(3).EdgeColor = 'none';
% 
% hyobs = slice(z,x,y,abs(Wo),[],[],yobs); 
% hyobs(1).FaceColor = 'interp'; 
% hyobs(2).FaceColor = 'interp';
% hyobs(3).FaceColor = 'interp';
% hyobs(1).EdgeColor = 'none'; 
% hyobs(2).EdgeColor = 'none'; 
% hyobs(3).EdgeColor = 'none'; 


daspect([1.5,.1,.1]) 
axis tight 
view(44,16) 
camzoom(1.4) 
camproj perspective
axis off



%distances of propagation
zprop  = [0, 6*Dz/3 , Dz/2, 0.99*Dz];
hzprop = slice(znorm,xnorm,ynorm,abs(bb),zprop,[],[]); 
hzprop(1).FaceColor = 'interp'; 
hzprop(2).FaceColor = 'interp';
hzprop(3).FaceColor = 'interp';
hzprop(4).FaceColor = 'interp';
hzprop(1).EdgeColor = 'none'; 
hzprop(2).EdgeColor = 'none'; 
hzprop(3).EdgeColor = 'none'; 
hzprop(4).EdgeColor = 'none'; 
colormap(mapgreen)


% hzprop(3).AlphaData = abs(Wo(:,64,:));
% hzprop(4).AlphaData = abs(Wo(:,128,:));
hold off
%% values of Alpha
% planes rotates 45 -45
% h45.FaceAlpha       = 0.0;
% h_45.FaceAlpha      = 0.0;
% % slices on borders 
% hxmin.FaceAlpha     = 1; 
% hxmax.FaceAlpha     = 1; 
% hy.FaceAlpha        = 0;
% hz.FaceAlpha        = 0; 
% slices on obstruction
% hxobs(1).FaceAlpha  = 0.0;
% hxobs(2).FaceAlpha  = 0.0; % x axis
% hxobs(3).FaceAlpha  = 0.0;
% hyobs(1).FaceAlpha  = 0.0;
% hyobs(2).FaceAlpha  = 0.0;
% hyobs(3).FaceAlpha  = 0.0;

hzprop(1).FaceAlpha = 0.8;
hzprop(2).FaceAlpha = 0;
hzprop(3).FaceAlpha = 0.8;
hzprop(4).FaceAlpha = 0.8;

xlim([0 1.1*Dz])
 set(gcf, 'color', [0 0 0])
set(gca,'XColor','g')
set(gca,'YColor','g')
set(gca,'ZColor','g')
set(gca, 'color', [0 0 0])
axis on
xlabel('$\zeta$')
ylabel('$\eta$')
zlabel('$\xi$')
% set(gcf,'Color',[1 1 1])
%% w
fieldfw = (W(:,floor(Nz),:));
fieldfw = reshape(fieldfw,[Nx,Nx]);
figure(2)
plotOpticalField(x,x,abs(fieldfw),mapgreen,'n')
axis off

hold on
plotPropagatedRaysInZ(rayH11,rayH12,rayH21,rayH22,1,1,floor(Nz))
  hold off
%%
%% w
fieldf = (Wo(:,floor(Nz),:));
fieldf = reshape(fieldf,[Nx,Nx]);
figure(2)
plotOpticalField(x,x,abs(fieldf),mapgreen,'n')
axis off

hold on
plotPropagatedRaysInZ(rayH11,rayH12,rayH21,rayH22,1,1,floor(Nz))
  hold off
  

  %%
  g  = fieldf(512,:);
  gw = fieldfw(512,:);
  figure(3)
  plot(x/InitialWaist,abs(g),'.')
  hold on
  plot(x/InitialWaist,abs(gw),'-')
  hold off
  legend('Obstructed Hermite Gauss','Hermite Gauss')
  title('z = (Rayleigh Distance)/2')
  xlabel('x at y = 0.5*\sigma_y')
  ylabel('Amplitude')
  
  %%
  dtx = 00;
  
  g  = fieldf(512+dtx,:);
  gw = fieldfw(512+dtx,:);
%   close(figure(3))
  figure(3)
  p(1) = plot(x,abs(g),'.');
  hold on
  p(2) = plot(x,abs(gw),'-');
  hold off
  set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
  set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

  title('z = (Rayleigh Distance)/2')
  xlabel('x at y = 0')
  ylabel('Amplitude')
  xline(rayH11(end).xCoordinate(1))
  xline(-rayH11(end).xCoordinate(1))
  xline(rayH22(end).xCoordinate(1))
  xline(-rayH22(end).xCoordinate(1))
    legend([p(1),p(2)],{'Obstructed Hermite Gauss','Hermite Gauss'})
%   set(get(get(p(3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%%
  dtx = 0;
  
  g  = fieldf(:,512+dtx);
  gw = fieldfw(:,512+dtx);
%   close(figure(3))
  figure(3)
  p(1) = plot(x,abs(g),'.');
  hold on
  p(2) = plot(x,abs(gw),'-');
  hold off
  set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
  set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

  title('z = (Rayleigh Distance)/2')
  xlabel('y at x = 0')
  ylabel('Amplitude')
  xline(rayH11(end).yCoordinate(3))
  xline(-rayH11(end).yCoordinate(3))
  xline(rayH22(end).yCoordinate(3))
  xline(-rayH22(end).yCoordinate(3))
    legend([p(1),p(2)],{'Obstructed Hermite Gauss','Hermite Gauss'})
    
    
%% last figure on paper
 
morado = [0.3010, 0.7450, 0.9330];

close(figure(2))
fig2 = figure(2); 
% set(gca,'Color','k')
set(gcf, 'color', [0 0 0])

fig2.Position = [537 381 1106 597];
ha  = tight_subplot(2,2,[.01 .01],[.1 .1],[.1 .07]);

dtx1 = 135;

fieldf  = (Wo(:,floor(Nz),:));
fieldf  = reshape(fieldf,[Nx,Nx]);

fieldfw = (W(:,floor(Nz),:));
fieldfw = reshape(fieldfw,[Nx,Nx]);


axes(ha(1))
set(gca,'Color','k')
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(fieldf),mapgreen,'n')

% plotPropagatedRaysInZ(rayH11,rayH12,rayH21,rayH22,1,1,floor(Nz))
% xline(x(512+dtx1),':','color','r','LineWidth',2)

hold on
plotPropagatedRaysInZ(rayH11,rayH12,rayH21,rayH22,InitialWaist,InitialWaist,z_index)
hold on
dtx0 = -20;
g0   = fieldf(512+dtx0,:);
p(1) = plot(x/InitialWaist,80*abs(g0)+dtx0/InitialWaist,'.','color','r');
yline(dtx0/InitialWaist,':','color','r','LineWidth',2)
ha(1).Color = 'blue';

g1   = fieldf(512+dtx1,:);
p(1) = plot(x/InitialWaist,100*abs(g1)+dtx1/InitialWaist,'.','color','b');
yline(dtx1/InitialWaist,':','color','b','LineWidth',2)
% p(2) = plot(x,abs(gw),'-');

hold off


axes(ha(3))
set(gca,'Color','k')
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(fieldfw),mapgreen,'n')

hold on
% plotPropagatedRaysInZ(rayH11,rayH12,rayH21,rayH22,InitialWaist,InitialWaist,z_index)
hold on

gw0   = fieldfw(512+dtx0,:);
p(1) = plot(x/InitialWaist,100*abs(gw0)+dtx0/InitialWaist,'.','color',morado);
yline(dtx0/InitialWaist,':','color','b','LineWidth',2)


gw1   = fieldfw(512+dtx1,:);
p(1) = plot(x/InitialWaist,100*abs(gw1)+dtx1/InitialWaist,'.','color','m');
yline(dtx1/InitialWaist,':','color','m','LineWidth',2)

hold off

axes(ha(4))

plot(x/InitialWaist,abs(gw0),'color',morado,'LineWidth',2)
hold on
plot(x/InitialWaist,abs(g0),'-.','color','r','LineWidth',2)
pbaspect([2.5 1 1])
% xline( rayH11(end).xCoordinate(1)/InitialWaist,':','LineWidth',1.7,'Color','m')
% xline(-rayH11(end).xCoordinate(1)/InitialWaist,':','LineWidth',1.7,'Color','c')
% xline( rayH22(end).xCoordinate(1)/InitialWaist,':','LineWidth',1.7,'Color','c')
% xline(-rayH22(end).xCoordinate(1)/InitialWaist,':','LineWidth',1.7,'Color','m')
hold off
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('Amplitude')
set(gca,'Color','k')
set(gcf,'color', [0 0 0])

axes(ha(2))
plot(x/InitialWaist,abs(gw1),'color','m','LineWidth',2)
hold on
plot(x/InitialWaist,abs(g1),'-.','color','b','LineWidth',2)
xline( rayH11(end).xCoordinate(1)/InitialWaist,':','LineWidth',1.7,'Color','m')
xline(-rayH11(end).xCoordinate(1)/InitialWaist,':','LineWidth',1.7,'Color','c')
xline( rayH22(end).xCoordinate(1)/InitialWaist,':','LineWidth',1.7,'Color','c')
xline(-rayH22(end).xCoordinate(1)/InitialWaist,':','LineWidth',1.7,'Color','m')
pbaspect([2.5 1 1])
hold off
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('Amplitude')
set(gca,'Color','k')
set(gcf,'color', [0 0 0])

% set(ha(1),'ycolor','w')
% set(ha(2),'ycolor','w')
% set(ha(3),'ycolor','w')
% set(ha(4),'ycolor','w')
% set(ha(1),'xcolor','w')
% set(ha(2),'xcolor','w')
% set(ha(3),'xcolor','w')
% set(ha(4),'xcolor','w')

export_fig('Selfhealingv5','-png','-transparent')

%%
dtx  = 00;
g    = fieldf(512+dtx,:);
gw   = fieldfw(512+dtx,:);
p(1) = plot(x,abs(g),'.');
hold on
p(2) = plot(x,abs(gw),'-');
hold off
set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%title('$z = \frac{z_R}{2}$','Interpreter','latex','FontSize',18)
xlabel('$y$ at $x = 0$','Interpreter','latex','FontSize',18)
ylabel('Amplitude')
xline(rayH11(end).xCoordinate(1),':','color','r','LineWidth',2)
xline(-rayH11(end).xCoordinate(1),':','color','m','LineWidth',2)
xline(rayH22(end).xCoordinate(1),':','color','m','LineWidth',2)
xline(-rayH22(end).xCoordinate(1),':','color','r','LineWidth',2)
%legend([p(1),p(2)],{'Obstructed Hermite Gauss','Hermite Gauss'})

axes(ha(3))
g    = fieldf(512+dtx1,:);
gw   = fieldfw(512+dtx1,:);
p(1) = plot(x,abs(g),'.');
hold on
p(2) = plot(x,abs(gw),'-');

set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%title('$z = \frac{z_R}{2}$','Interpreter','latex','FontSize',18)
xlabel(['$y$ at $x =',num2str(x(512+dtx1)),'$'],'Interpreter','latex','FontSize',18)
ylabel('Amplitude')
xline(rayH11(end).xCoordinate(1),':','color','r','LineWidth',2)
xline(-rayH11(end).xCoordinate(1),':','color','m','LineWidth',2)
xline(rayH22(end).xCoordinate(1),':','color','m','LineWidth',2)
xline(-rayH22(end).xCoordinate(1),':','color','r','LineWidth',2)
legend([p(1),p(2)],{'Obstructed Hermite Gauss','Hermite Gauss'})
hold off