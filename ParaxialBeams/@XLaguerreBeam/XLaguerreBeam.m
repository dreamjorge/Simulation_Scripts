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
    waistL  = waistLaguerre(zDistance,InitialWaist,RayleighDistance,nu,mu);
  end
  
  methods
    
    function PhiPhase = get.PhiPhase(obj)
      PhiPhase = (obj.l+obj.p).*obj.GouyPhase;
    end

    function LaguerreWaist = get.LaguerreWaist(obj)
      LaguerreWaist = sqrt(obj.l+obj.p);
    end
    
    function LaguerreAmplitude = get.LaguerreAmplitude(obj)
      LaguerreAmplitude = (sqrt(2)*(obj.RadialCoordinate)./obj.Waist).^(obj.l);
    end

    function Normalization = get.Normalization(obj)
      Normalization =  sqrt(2*factorial(obj.p)/(pi*factorial(obj.p+abs(obj.l))));
    end

    
    function Laguerre = XLaguerreBeam(x,y,PropagationDistance, InitialWaist,Wavelength,l,p)

      if nargin == 7    
        

        super_args{1} = x;
        super_args{2} = y;
        super_args{3} = PropagationDistance;
        super_args{4} = InitialWaist;
        super_args{5} = Wavelength;
      else
        
      end
      
      %Laguerre@GaussianParameters(super_args{3:end});
      Laguerre@GaussianBeam(super_args{:}); 
      
      Laguerre.l = l;
      Laguerre.p = p;
      [Laguerre.ThetaCoordinate,Laguerre.RadialCoordinate] = cart2pol(x,y);
      
      %% Optical Field
      Laguerre.OpticalField = Laguerre.Normalization.*...
                              Laguerre.LaguerreAmplitude.*... 
                              exp(-1i*Laguerre.PhiPhase).*exp(1i*abs(l)*Laguerre.ThetaCoordinate).*...
                              XLaguerreBeam.XAssociatedLaguerrePolynomial(l,abs(p),(2*Laguerre.RadialCoordinate.^2)./Laguerre.Waist.^2).*...
                              Laguerre.OpticalField;
    end
  end
  
  
  
  
  
end