function [XLGB]=xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,r,th,z) 
    
    k    = 2*zo/wo^2;

    wz   = wo*sqrt(1+(z./zo).^2);

    thr  = 45;
    
    XLGB = (2*r.^2./(wz.^2)).^(mu/2).*XAssociatedLaguerrePolynomial(thr,nu,mu,2*r.^2./(wz.^2)).*physicalRadialGaussianBeam(wo,zo,r,z).*exp(1i*mu*th);
    
end