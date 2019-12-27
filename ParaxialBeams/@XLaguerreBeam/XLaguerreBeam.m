classdef XLaguerreBeam < GaussianBeam

  properties
    l
    p
  end
  
  properties (Dependent)
    LaguerreWaist
    PhiPhase
    LaguerreAmplitude
    Normalization
  end
  
  properties (Hidden)
    RadialCoordinate
    ThetaCoordinate
  end
  
  methods(Static)
    XLg     = XAssociatedLaguerrePolynomial(nu,mu,x);
  end
  
  methods
    
    function PhiPhase = get.PhiPhase(obj)
      PhiPhase = (abs(obj.l)+2*obj.p).*obj.GouyPhase;
    end

    function LaguerreWaist = get.LaguerreWaist(obj)
      LaguerreWaist = sqrt(obj.l+obj.p);
    end
    
    function LaguerreAmplitude = get.LaguerreAmplitude(obj)
      LaguerreAmplitude = 1;...(sqrt(2)*(obj.RadialCoordinate)./obj.Waist).^(obj.l);
    end

    function Normalization = get.Normalization(obj)
      Normalization =  sqrt(2*factorial(obj.p)/(pi*factorial(obj.p+abs(obj.l))));
    end

    
    function XLaguerre = XLaguerreBeam(x,y,PropagationDistance, InitialWaist,Wavelength,nu,mu)

      if nargin == 7    
        super_args{1} = x;
        super_args{2} = y;
        super_args{3} = PropagationDistance;
        super_args{4} = InitialWaist;
        super_args{5} = Wavelength;
      else
        
      end
      
      %Laguerre@GaussianParameters(super_args{3:end});
      XLaguerre@GaussianBeam(super_args{:}); 
      
      XLaguerre.l = nu;
      XLaguerre.p = mu;
      [XLaguerre.ThetaCoordinate,XLaguerre.RadialCoordinate] = cart2pol(x,y);
      
      %% Optical Field
      XLaguerre.OpticalField =... XLaguerre.Normalization.*...
                              ...XLaguerre.LaguerreAmplitude.*... 
                              exp(1i*XLaguerre.PhiPhase).*exp(-1i*abs(nu)*XLaguerre.ThetaCoordinate).*...
                              XLaguerreBeam.XAssociatedLaguerrePolynomial(nu,abs(mu),(2*XLaguerre.RadialCoordinate.^2)./XLaguerre.Waist.^2).*...
                              XLaguerre.OpticalField;
    end
  end
  
  
  
  
  
end