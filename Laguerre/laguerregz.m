function [Ln]=laguerregz(nu,mu,wo,zl,x,z) 
    k=2*zl/wo^2;
    Rz=z+zl^2./z;
    wz=wo*sqrt(1+(z./zl).^2);
    phiz=(2*nu+mu+1)*atan(z./zl);
%     Ao=sqrt(2*factorial(nu)/(1+kroneckerDelta(mu))*pi*factorial(nu+mu));
Ln=LaguerreG(nu,abs(mu),2*x.^2./(wz.^2)).*exp(1i*(k*x.^2./(2*Rz)-phiz));%./(wz.^(mu+1));
end