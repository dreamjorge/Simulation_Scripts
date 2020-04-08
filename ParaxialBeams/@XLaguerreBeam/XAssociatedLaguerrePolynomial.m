function [Xnmr]= XAssociatedLaguerrePolynomial(n,m,x)
%This function calculates second solution of associated laguerre equation
%for n,m integers.
nrterms = 48;

  % first term
  a1 = x.^(n+1)/(gamma(n+2)*gamma(m+1+n+1));

  for k = n+2:nrterms % summing terms

      a1 = a1+gamma(k-n)*x.^k/(gamma(k+1)*gamma(m+1+k));

  end
  % -a1, this term is added to solution.

  a2 = x.^(-1)/(gamma(n+2)*gamma(m));

  for k = 2:m % summing terms

      a2 = a2+gamma(k)*x.^(-k)/(gamma(n+k+1)*gamma(m-k+1));

  end

  % (-1)^n*a2

  a3 = (-log(x)+psi(1)+psi(m+1)-psi(n+1))...
           /(gamma(n+1)*gamma(m+1));

  for k = 1:n % summing terms

      a3=a3+(-log(x)+psi(k+1)+psi(m+k+1)-psi(n-k+1))...
           .*(-x).^k/(gamma(n-k+1)*gamma(m+k+1)*gamma(k+1));

  end
  % (-1)^(n)*a3

  % Second Solution
  Xnmr = (-1)^(n+1)/((n+(m+1)/2)^(m/2))*gamma(n+m+1)/pi*(-a1+((-1)^n)*(a2+a3));

end