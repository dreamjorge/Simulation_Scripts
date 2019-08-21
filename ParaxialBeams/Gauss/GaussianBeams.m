% color of beam
mapgreen = AdvancedColormap('kg',256,[0  255]/255);  %color del haz

%-----------------------Normalized Gaussian Beams------------------------%

Ns = 2^9;                     % number of points
ns = -Ns/2:Ns/2-1;            % index vector with 1 of resolution

%%vector s
Ds = 6;                       % size of window of vector
ds = Ds/Ns;                   % resolution of vector
s  = ns.*ds;                  % vector

N = 2^9;                      % number of points
n = -N/2:N/2-1;               % index vector with 1 of resolution

%%vector xi
Dxi = 5;                      % size of window
dxi = Dxi/N;                  % resolution of vector
xi  = n.*dxi;                 % vector

%%vector eta
Deta = 5;                     % size of window
deta = Deta/N;                % resolution of vector
eta  = n.*deta;               % vector

[Xi,Eta] = meshgrid(xi,eta);

%%Parameters of gaussian beam

wo     = 1;
sigmao = 2*wo;
so     = wo^2;    



u0   = 1;
Rs   = s./2 + 2./s;          % Radius of curvature of beam
ws   = sqrt(1+s.^2/4);       % Waist of beam
Phis = atan(s);              % Phase of beam

%% lateral beam

uxi  = zeros(N,Ns);
ueta = zeros(N,Ns);
us  = zeros(N,N,Ns);

for ii = 1:Ns %% run all points in s 
    
    % Calculating lateral beam
    ueta(:,ii) = (u0./ws(ii)).*exp(1i*(xi(N/2+1).^2+eta.^2)/(2*Rs(ii))).*exp(-(xi(N/2+1).^2+eta.^2)/(2*(ws(ii).^2))).*exp(-1i*Phis(ii));
    uxi(:,ii)  = (u0./ws(ii)).*exp(1i*(xi.^2+eta(N/2+1).^2)/(2*Rs(ii))).*exp(-(xi.^2+eta(N/2+1).^2)/(2*(ws(ii).^2))).*exp(-1i*Phis(ii)); 
    us(:,:,ii)  = (u0./ws(ii)).*exp(1i*(Xi.^2+Eta.^2)./(2*Rs(ii))).*exp(-(Xi.^2+Eta.^2)/(2*(ws(ii).^2))).*exp(-1i*Phis(ii)); 
end
%% in s = 0

indexS = find(s==0);

uso   = (u0./ws(indexS)).*exp(1i*(Xi.^2+Eta.^2)./(2*Rs(indexS))).*exp(-(Xi.^2+Eta.^2)/(2*(ws(indexS).^2))).*exp(-1i*Phis(indexS)); 

usoxi =  (u0./ws(indexS)).*exp(1i*(xi.^2+eta(N/2+1).^2)/(2*Rs(indexS))).*exp(-(xi.^2+eta(N/2+1).^2)/(2*(ws(indexS).^2))).*exp(-1i*Phis(indexS)); 


sigma = ws(indexS);
sigmay =  (u0./ws(indexS)).*exp(1i*(ws(indexS).^2+eta(N/2+1).^2)/(2*Rs(indexS))).*exp(-(ws(indexS).^2+eta(N/2+1).^2)/(2*(ws(indexS).^2))).*exp(-1i*Phis(indexS)); 
sigmay2 =  (u0./(4*ws(indexS))).*exp(1i*((2*ws(indexS)).^2+eta(N/2+1).^2)/(2*Rs(indexS))).*exp(-((2*ws(indexS)).^2+eta(N/2+1).^2)/(2*((2*ws(indexS)).^2))).*exp(-1i*Phis(indexS)); 


%%daspect([100 1 1])


figure(1)

subplot(1,2,1)
imagesc(xi,eta,abs(uso))
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
set(gca,'YDir','normal')
colormap(mapgreen)
daspect([1 1 1])
xticksv =[-2*sigma -sigma  0 sigma 2*sigma];
xticklabelsv={'$-2\sigma(z)$','$-\sigma(z)$','$0$','$\sigma(z)$','$2\sigma(z)$'};
set(gca,'xtick',xticksv); 
set(gca,'xticklabel',xticklabelsv)
set(gca,'ytick',xticksv); 
set(gca,'yticklabel',xticklabelsv)




set(groot,'defaultAxesTickLabelInterpreter','latex');  
subplot(1,2,2)
plot(xi,abs(usoxi).^2)

xlabel('$x$','Interpreter','latex')
ylabel('Amplitude')
set(gca,'YTickLabel',[]);
%%set(gca,'XTickLabel',[]);
xticksv =[-sqrt(2)*2*sigma  -2*sigma -sigma  0 sigma 2*sigma sqrt(2)*2*sigma];
xticklabelsv={'$-2\sqrt{2}\sigma(z)$','$-2\sigma(z)$','$-\sigma(z)$','$0$','$\sigma(z)$','$2\sigma(z)$','$2\sqrt{2}\sigma(z)$'};
set(gca,'xtick',xticksv); 
set(gca,'xticklabel',xticklabelsv);
p1 = line([sigma sigma], [0 abs(sigmay).^2],'Color','r');
p2 = line([-sigma -sigma], [0 abs(sigmay).^2],'Color','r');
p3 = line([-2*sigma -2*sigma], [0 abs(sigmay2).^2],'Color','b');
p4 = line([2*sigma 2*sigma], [0 abs(sigmay2).^2],'Color','b');
p5 = line([sqrt(2)*2*sigma sqrt(2)*2*sigma], [0 abs(sigmay2).^2],'Color','b');
p6 = line([sqrt(2)*2*sigma sqrt(2)*2*sigma], [0 abs(sigmay2).^2],'Color','b');
suptitle('Gaussian Intensity Profile of Gaussian Beam')



%%
uso  = us(:,:,1);
usxo = us(:,N/2+1,1);


%%

figure(2)
imagesc(s,xi,abs(uxi).^2)
set(gca,'YDir','normal')
colormap(mapgreen)
hold on
p1 = plot(s ,ws,'r','LineWidth',1.5);
p2 = plot(s,-ws,'r','LineWidth',1.5);
p3 = plot(s, Rs,'LineWidth',1.5)
p4 = plot(s, atan(1/2)*s,'c','LineWidth',1.5);
hold off
xlabel('$s [z/z_R]$','Interpreter','latex')
ylabel('$\eta [x/w_0]$','Interpreter','latex')
title('Lateral Beam')
axis square
daspect([1 1 1])
legend([p1 p3 p4],{'$w(s)$','$R(s)$','$w(s)=\frac{1}{2}s$'},'Interpreter','latex')


figure(3)
imagesc(s,xi,angle(uxi))

%%

% us = (u0./ws).*exp(1i*(xi.^2+eta.^2)/(2*Rs)).*exp(-(xi.^2+eta.^2)/(2*(ws.^2))).*exp(-1i*Phis);
% 
% figure(1)
% plot(s)
% 
% 
% 
% figure(2)
% plot(ws)