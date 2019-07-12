function [Xn]= XLaguerreG(r,n,m,x)

Xn=(-1)^(n+1)/(gamma(n+1)*gamma(m+1)*(n+(m+1)/2)^(m/2))...
       *exp(-x/2).*x.^(m/2).*XLaguerre(r,n,m,x);

end