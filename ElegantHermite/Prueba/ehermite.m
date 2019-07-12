function [HG,NHG]= ehermite(nu,x)
%Funcion que determina las series de Hermite para un eigenvalor nu
%fpar da la serie para el eigenvalor nu
%fimpar da la serie para el eigenvalor nu

%hacemos a_0=1, y a_1=1 del desarrollo en serie para simplificar
an=1;
bn=1;

%Los dos primeros terminos de las series, el del termino par e impar
fpar=1;
fimpar=x;

%Determinando los terminos de la serie hasta un termino n, a mayor n mayor
%precisión
n=floor(nu+nu/2);
for k = 0:n
 %Serie par   
 an=an*(2*((2*k)-nu))/(((2*k)+1)*((2*k)+2));        %termino a_{2n} par de la serie, si n=nu/2 el termino a_{2n} se hace cero 
 fpar=fpar+an*(x).^(2*k+2);                         %realizando la serie hasta el termino nu, pues cuando n=nu/2 los demas terminos son cero 
 %Serie impar
 bn=bn*(2*((2*k+1)-nu))/(((2*k+1)+1)*((2*k+1)+2));  %termino a_{2n+1} impar de la serie, si n=(nu-1)/2 el termino a_{2n+1} se hace cero 
 fimpar=fimpar+bn*(x).^(2*k+3);                     %realizando la serie hasta el termino nu, pues cuando  n=(nu-1)/2 los demas terminos son cero 
end

%factor (-1) en la serie, el cual intercambia el signo de la misma en el
%eje y
% fpar=(1/sqrt((2^nu)*factorial(nu)*sqrt(pi)))*(-1)^(nu/2)*fpar;
% fimpar=(1/sqrt((2^nu)*factorial(nu)*sqrt(pi)))*(-1)^((nu-1)/2)*fimpar;

% fpar=(-1)^(nu/2)*fpar;
% fimpar=(-1)^((nu-1)/2)*fimpar;

%multiplicando por el factor gaussiano para obtener las funciones
%Hermite Gauss

fpar=fpar.*exp(-(x.^2));
fimpar=fimpar.*exp(-(x.^2));

Norma=sqrt(2*nu+1);

if(mod(nu,2)~=0)        %Condicional para nu impar

% yfpar=imag(fpar);
fimpar=Norma*fimpar;

HG=fimpar;
NHG=fpar;
else                    %Condicional cuando nu no es impar
fpar=fpar;
fimpar=Norma*(fimpar);
% fimpar=Norma*imag(fimpar);   

NHG=fimpar;
HG=fpar;
end
end