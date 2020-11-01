 ray11 = struct();
  ray12 = struct();
  ray21 = struct();
  ray22 = struct();
  for z_index = 1:Nz

        for ray_index = 1:TotalRays
          
          ray11(ray_index).z(z_index) = rayH11(z_index).zCoordinate(ray_index);
          ray11(ray_index).x(z_index) = rayH11(z_index).xCoordinate(ray_index);
          ray11(ray_index).y(z_index) = rayH11(z_index).yCoordinate(ray_index);
          
          ray12(ray_index).z(z_index) = rayH12(z_index).zCoordinate(ray_index);
          ray12(ray_index).x(z_index) = rayH12(z_index).xCoordinate(ray_index);
          ray12(ray_index).y(z_index) = rayH12(z_index).yCoordinate(ray_index);
          
          ray21(ray_index).z(z_index) = rayH21(z_index).zCoordinate(ray_index);
          ray21(ray_index).x(z_index) = rayH21(z_index).xCoordinate(ray_index);
          ray21(ray_index).y(z_index) = rayH21(z_index).yCoordinate(ray_index);
          
          ray22(ray_index).z(z_index) = rayH22(z_index).zCoordinate(ray_index);
          ray22(ray_index).x(z_index) = rayH22(z_index).xCoordinate(ray_index);
          ray22(ray_index).y(z_index) = rayH22(z_index).yCoordinate(ray_index);
          
        end
 
  end    

%%

gxt = Wo(512,:,:);
  
  gtx = squeeze(gxt);
  
  fig6 = figure(6);
  fig6.Position = [382 228 1375 537];
  pcolor(z/RayleighDistance,x/InitialWaist, abs(gtx').^2);
  shading interp
  xlabel('$z/z_R$','Interpreter','latex','FontSize',18)
  ylabel('$y/w_o$','Interpreter','latex','FontSize',18)
  xlim([0 0.5])
  export_fig('HermiteLateralX','-png','-transparent')
  hold on
  indexray = 4;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
  indexray = 2;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
   
%    plot(z/RayleighDistance,-HPz.HermiteWaistX/(2*InitialWaist),'-.','Linewidth',2,'Color','r')
  hold off
  export_fig('HermiteLateralXRays','-png','-transparent')
  %%
  gxt = Wo(:,:,512);
  
  gtx = squeeze(gxt);
  
  fig6 = figure(6);
  fig6.Position = [382 228 1375 537];
  pcolor(z/RayleighDistance,x/InitialWaist, abs(gtx).^2);
  shading interp
  xlabel('$z/z_R$','Interpreter','latex','FontSize',18)
  ylabel('$y/w_o$','Interpreter','latex','FontSize',18)
  xlim([0 0.5])
  export_fig('HermiteLateralY','-png','-transparent')
  hold on
  indexray = 4;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
  indexray = 2;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
   
%    plot(z/RayleighDistance,-HPz.HermiteWaistX/(2*InitialWaist),'-.','Linewidth',2,'Color','r')
  hold off
  export_fig('HermiteLateralYRays','-png','-transparent')

  
  
  
  
  
 %% 
  gxt = W11o(512,:,:);
  
  gtx = squeeze(gxt);
  
  fig6 = figure(6);
  fig6.Position = [382 228 1375 537];
  pcolor(z/RayleighDistance,x/InitialWaist, abs(gtx').^2);
  shading interp
  xlabel('$z/z_R$','Interpreter','latex','FontSize',18)
  ylabel('$y/w_o$','Interpreter','latex','FontSize',18)
  xlim([0 0.5])
  export_fig('HermiteLateralX11','-png','-transparent')
  hold on
  indexray = 4;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
  indexray = 2;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
   
%    plot(z/RayleighDistance,-HPz.HermiteWaistX/(2*InitialWaist),'-.','Linewidth',2,'Color','r')
  hold off
  export_fig('HermiteLateralX11Rays','-png','-transparent')
  
  %%
    gxt = W11o(:,:,512);
  
  gtx = squeeze(gxt);
  
  fig6 = figure(6);
  fig6.Position = [382 228 1375 537];
  pcolor(z/RayleighDistance,x/InitialWaist, abs(gtx).^2);
  shading interp
  xlabel('$z/z_R$','Interpreter','latex','FontSize',18)
  ylabel('$y/w_o$','Interpreter','latex','FontSize',18)
  xlim([0 0.5])
  export_fig('HermiteLateralY11','-png','-transparent')
  hold on
  indexray = 4;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
  indexray = 2;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
   
%    plot(z/RayleighDistance,-HPz.HermiteWaistX/(2*InitialWaist),'-.','Linewidth',2,'Color','r')
  hold off
  export_fig('HermiteLateralY11Rays','-png','-transparent')
  
  
  
%%

  gxt = W22o(512,:,:);
  
  gtx = squeeze(gxt);
  
  fig6 = figure(6);
  fig6.Position = [382 228 1375 537];
  pcolor(z/RayleighDistance,x/InitialWaist, abs(gtx').^2);
  shading interp
  xlabel('$z/z_R$','Interpreter','latex','FontSize',18)
  ylabel('$y/w_o$','Interpreter','latex','FontSize',18)
  xlim([0 0.5])
  export_fig('HermiteLateralX22','-png','-transparent')
  hold on
  indexray = 4;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
  indexray = 2;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
   
%    plot(z/RayleighDistance,-HPz.HermiteWaistX/(2*InitialWaist),'-.','Linewidth',2,'Color','r')
  hold off
  export_fig('HermiteLateralX22Rays','-png','-transparent')
  
  %%
   gxt = W22o(:,:,512);
  
  gtx = squeeze(gxt);
  
  fig6 = figure(6);
  fig6.Position = [382 228 1375 537];
  pcolor(z/RayleighDistance,x/InitialWaist, abs(gtx).^2);
  shading interp
  xlabel('$z/z_R$','Interpreter','latex','FontSize',18)
  ylabel('$y/w_o$','Interpreter','latex','FontSize',18)
  xlim([0 0.5])
  export_fig('HermiteLateralY22','-png','-transparent')
  hold on
  indexray = 4;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
  indexray = 2;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
   
%    plot(z/RayleighDistance,-HPz.HermiteWaistX/(2*InitialWaist),'-.','Linewidth',2,'Color','r')
  hold off
  export_fig('HermiteLateralY22Rays','-png','-transparent')
  
  
  
  %%
  
  gxt = W12o(512,:,:);
  
  gtx = squeeze(gxt);
  
  fig6 = figure(6);
  fig6.Position = [382 228 1375 537];
  pcolor(z/RayleighDistance,x/InitialWaist, abs(gtx').^2);
  shading interp
  xlabel('$z/z_R$','Interpreter','latex','FontSize',18)
  ylabel('$y/w_o$','Interpreter','latex','FontSize',18)
  export_fig('HermiteLateralX12','-png','-transparent')
  hold on
  indexray = 4;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
  indexray = 2;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
   
%    plot(z/RayleighDistance,-HPz.HermiteWaistX/(2*InitialWaist),'-.','Linewidth',2,'Color','r')
  hold off
  export_fig('HermiteLateralX12Rays','-png','-transparent')
  
 %%
   gxt = W12o(:,:,512);
  
  gtx = squeeze(gxt);
  
  fig6 = figure(6);
  fig6.Position = [382 228 1375 537];
  pcolor(z/RayleighDistance,x/InitialWaist, abs(gtx).^2);
  shading interp
  xlabel('$z/z_R$','Interpreter','latex','FontSize',18)
  ylabel('$y/w_o$','Interpreter','latex','FontSize',18)
  xlim([0 0.5])
  export_fig('HermiteLateralY12','-png','-transparent')
  hold on
  indexray = 4;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
  indexray = 2;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
   
%    plot(z/RayleighDistance,-HPz.HermiteWaistX/(2*InitialWaist),'-.','Linewidth',2,'Color','r')
  hold off
  export_fig('HermiteLateralY12Rays','-png','-transparent')
  
  
 %%
    
  gxt = W21o(512,:,:);
  
  gtx = squeeze(gxt);
  
  fig6 = figure(6);
  fig6.Position = [382 228 1375 537];
  pcolor(z/RayleighDistance,x/InitialWaist, abs(gtx').^2);
  shading interp
  xlabel('$z/z_R$','Interpreter','latex','FontSize',18)
  ylabel('$y/w_o$','Interpreter','latex','FontSize',18)
  xlim([0 0.5])
  export_fig('HermiteLateralX21','-png','-transparent')
  hold on
  indexray = 4;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
  indexray = 2;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
   
%    plot(z/RayleighDistance,-HPz.HermiteWaistX/(2*InitialWaist),'-.','Linewidth',2,'Color','r')
  hold off
  export_fig('HermiteLateralX21Rays','-png','-transparent')
  
  
  
  
   %%
    
  gxt = W21o(:,:,512);
  
  gtx = squeeze(gxt);
  
  fig6 = figure(6);
  fig6.Position = [382 228 1375 537];
  pcolor(z/RayleighDistance,x/InitialWaist, abs(gtx).^2);
  shading interp
  xlabel('$z/z_R$','Interpreter','latex','FontSize',18)
  ylabel('$y/w_o$','Interpreter','latex','FontSize',18)
  xlim([0 0.5])
  export_fig('HermiteLateralY21','-png','-transparent')
  hold on
  indexray = 4;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
  indexray = 2;
   plot(ray11(indexray).z/RayleighDistance, ray11(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','r')
   plot(ray22(indexray).z/RayleighDistance, ray22(indexray).x./InitialWaist,'-.','Linewidth',2,'Color','c')
   
%    plot(z/RayleighDistance,-HPz.HermiteWaistX/(2*InitialWaist),'-.','Linewidth',2,'Color','r')
  hold off
  export_fig('HermiteLateralY21Rays','-png','-transparent')