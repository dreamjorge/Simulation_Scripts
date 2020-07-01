function [Ln] = AssociatedLaguerrePolynomial(n,m,x)
%This function calculates Asociated Laguerre Polinomyal with recursive
%method (Alfonso formula)

    if (n == 0)

        Ln = 1;

    else

        if (n == 1)

            Ln = 1+m-x;

        else

            l0 = 1;
            l1 = 1+m-x;

            for k=2:n

                Ln = ((2*(k-1)+m+1-x).*l1-(k+m-1).*l0)./k;
                l0 = l1;
                l1 = Ln;

            end

        end

        
        Ln=(gamma(n+1)/gamma(n+m+1))./(((n+(m+1)/2)*x).^(-m/2)).*Ln;
        
    end

          
            
end 
    
