%%Hankel 11 slices
mapgreen = AdvancedColormap('kgg',256,[0 50 255]/255);  %color of beam


bb = W22o;
iindex = find(abs(bb)<=.000041);
 bb(iindex) = NaN;
close(figure(1))
fig = figure (1);
fig.Position = [680 365 772 613];

plotPropagatedRays3D(rayH12,1,1,1,'m');
plotPropagatedRaysSquareInZfor3D(rayH12,1,1,1,4,'m');
plotPropagatedRaysSquareInZfor3D(rayH12,1,1,1,17,'m');

hold on 

daspect([0.7,.1,.1]) 
axis tight 
view(44,16) 
camzoom(1.4) 
camproj perspective
axis off

zprop  = [0, (1/10)*RayleighDistance, (1/6)*RayleighDistance, (1/5)*RayleighDistance];
%zprop  = [0, 6*Dz/3 , Dz/2, 0.99*Dz];

hzprop = slice(z,x,y,abs(bb).^2,zprop,[],[]); 
hzprop(1).FaceColor = 'interp'; 
hzprop(2).FaceColor = 'interp';
hzprop(3).FaceColor = 'interp';
hzprop(4).FaceColor = 'interp';
hzprop(1).EdgeColor = 'none'; 
hzprop(2).EdgeColor = 'none'; 
hzprop(3).EdgeColor = 'none'; 
hzprop(4).EdgeColor = 'none'; 
colormap(mapgreen)

hold off
%
hzprop(1).FaceAlpha = 0.9;
hzprop(2).FaceAlpha = 0.9;
hzprop(3).FaceAlpha = 0.0;
hzprop(4).FaceAlpha = 0.9;

 xlim([0 1.02*(1/5)*RayleighDistance])
set(gcf, 'color', [0 0 0])
set(gca,'XColor','g')
set(gca,'YColor','g')
set(gca,'ZColor','g')
 set(gca, 'color', [0 0 0])
axis off
xlabel('\zeta')
ylabel('\eta')
zlabel('\xi')

export_fig('HankelHermite3D','-png','-transparent')

%% 

distances = 2*[25/12,10/2,inf];

texts     = {'6/25','1/10','0'};


FG = W22o;

close(figure(800))
fig800 = figure(800);
ha = tight_subplot(2,3,[.01 .01],[.05 .1],[.1 .07]);


kk = 1;
for jj = distances
  
  index = floor(Nz/jj);
  if index == 0
    index = 1;
  end
  

  gg1 = W22o(:,index,:);
  gg1 = reshape(gg1,[Nx,Nx]);
  figure(1000)
  plotOpticalField(x,x,abs(gg1).^2,mapgreen,'','');
  axis square
  axis off
  plotRays(rayH12(index),'m',1,1);
%   title(['$z$ = ', texts{kk}],'Interpreter','latex')

  kk = kk+1;

  export_fig(['HH11Propagation',num2str(jj)],'-png','-transparent')
%   saveas(gcf,['HH22Propagation',num2str(jj),'.png'])
end