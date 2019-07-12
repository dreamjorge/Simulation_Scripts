function [Xnmr]= XLaguerre(r,n,m,x)

%Empecemos por los términos que vamos a sumar
% La segunda solución metida en Mathematica está en Talachita de mi osita
% (carpeta Traveling Laguerre)
% % Usar vpa a partir de aproximadamente n=30
% x=vpa(0.01:.1:250,20);
% % x=1:10;
% n=vpa(50,20);
% % n=32;
% m=vpa(1,20);
% r=vpa(150,20);
% % for r=50:100
%%%%%%%%%%%%%%%%%%%%
% x=vpa(0.01:0.025:10,20);
% n=vpa(100,20);
% m=vpa(1,20);
% r=vpa(150,20);

%%%%%%%%%%%%%%%%%%%%

a1=x.^(n+1)/(gamma(n+2)*gamma(m+1+n+1));

for k=n+2:r        %Empezamos a calcular la suma
    
    a1=a1+gamma(k-n)*x.^k/(gamma(k+1)*gamma(m+1+k));
    
end
% -a1, el término que se debe sumar a la solución
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%

a2=x.^(-1)/(gamma(n+2)*gamma(m));

for k=2:m        %Empezamos a calcular la suma
    
    a2=a2+gamma(k)*x.^(-k)/(gamma(n+k+1)*gamma(m-k+1));
    
end
% (-1)^n*a2

%%%%%%%%%%%%%%%%%%%%

a3=(-log(x)+psi(1)+psi(m+1)-psi(n+1))...
         /(gamma(n+1)*gamma(m+1));

for k=1:n       %Empezamos a calcular la suma
    
    a3=a3+(-log(x)+psi(k+1)+psi(m+k+1)-psi(n-k+1))...
         .*(-x).^k/(gamma(n-k+1)*gamma(m+k+1)*gamma(k+1));
    
end
% (-1)^(n)*a3

Xnmr=gamma(n+1)*gamma(m+1)*gamma(n+m+1)/pi*(-a1+(-1)^n*(a2+a3));
end