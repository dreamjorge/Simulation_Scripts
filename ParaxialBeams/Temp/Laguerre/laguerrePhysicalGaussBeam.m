function [LGB] = laguerrePhysicalGaussBeam(nu,mu,wo,zo,r,th,z) 
%This function calculates Laguerre Gauss Beam with all parameters
    
    wz  = waistPhysicalGaussianBeam(z,wo,zo);

    LGB = AssociatedLaguerrePolynomial(nu,abs(mu),2*r.^2./(wz.^2)).*physicalRadialGaussianBeam(wo,zo,r,z).*exp(1i*mu*th);
    
end
