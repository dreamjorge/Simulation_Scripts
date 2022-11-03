function [HG,NHG] = getHermiteSolutions(nu,x)
% Solutions of Hermite Differential Equation in terms of Series

an=1;
bn=1;

fpar=1;
fimpar=x;

%value of n is for precision 
n=floor(nu+nu/2);
for k = 0:n
 %Serie par   
 an=an*(2*((2*k)-nu))/(((2*k)+1)*((2*k)+2));        
 fpar=fpar+an*(x).^(2*k+2);                        
 %Serie impar
 bn=bn*(2*((2*k+1)-nu))/(((2*k+1)+1)*((2*k+1)+2));  
 fimpar=fimpar+bn*(x).^(2*k+3);                     
end



Norma=sqrt(2*nu+1);

if(mod(nu,2)~=0)        

    % yfpar=imag(fpar);
    fimpar=Norma*fimpar;

    HG=fimpar;
    NHG=fpar;
else                   
    fpar=fpar;
    fimpar=Norma*(fimpar);


    NHG=fimpar;
    HG=fpar;
end
end