function [Xn]= XLaguerreG(r,n,m,x)

Xn=(gamma(n+1)/gamma(n+m+1))./(exp(x./2).*((n+(m+1)/2)*x).^(-m/2)).*...
                                 XLaguerre(r,n,m,x);

end