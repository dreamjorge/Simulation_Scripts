function [Xn]=xlaguerregz(nu,mu,wo,zl,x,z) 
    
k=2*zl/wo^2;
    Rz=z+zl^2./z;
    wz=wo*sqrt(1+(z./zl).^2);
    phiz=(2*nu+mu+1)*atan(z./zl);

Xn=XLaguerreG(45,nu,abs(mu),2*x.^2./(wz.^2)).*exp(1i*(k*x.^2./(2*Rz)-phiz));%./(wz.^(mu+1));
end