function [HG,NHG]= hermitegz(nu,wo,zl,x,z)
%Funcion que determina las series de Hermite para un eigenvalor nu
%fpar da la serie para el eigenvalor nu
%fimpar da la serie para el eigenvalor nu

% Parametros de Hermite
wz=wo.*sqrt(1+(z.^2)/(zl^2));
Rz=z+(zl^2./z);
Phizx=nu.*atan(z./zl);
A=wo./wz;
k=2*zl/wo^2;
%reescalando la funcion hermtite por wz
x=x./wz;


%hacemos a_0=1, y a_1=1 del desarrollo en serie para simplificar
an=1;
bn=1;

%Los dos primeros terminos de las series, el del termino par e impar
fpar=1;
fimpar=x;

%Determinando los terminos de la serie hasta un termino n, a mayor n mayor
%precisión
n=floor(nu+4*nu);
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

fpar=fpar.*exp(-(x.^2/2));
fimpar=fimpar.*exp(-(x.^2/2));

Norma=sqrt(2*nu+1);

if(mod(nu,2)~=0)        %Condicional para nu impar

fimpar=Norma*fimpar;

HG=A.*fimpar.*exp(1i*k*x.^2./(2*Rz)).*exp(1i*k*z).*exp(1i.*k.*Phizx);
NHG=A.*fpar.*exp(1i*k*x.^2./(2*Rz)).*exp(1i*k*z).*exp(1i.*k.*Phizx);
else                    %Condicional cuando nu no es impar
fimpar=Norma*(fimpar);
% fimpar=Norma*imag(fimpar);   

NHG=A.*fimpar.*exp(1i*k*x.^2./(2*Rz)).*exp(1i*k*z).*exp(1i.*k.*Phizx);
HG=A.*fpar.*exp(1i*k*x.^2./(2*Rz)).*exp(1i*k*z).*exp(1i.*k.*Phizx);


end




end