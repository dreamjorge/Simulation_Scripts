classdef ElegantLaguerreBeam
    % ElegantLaguerreBeam - Elegant Laguerre-Gaussian beam implementation
    % Consistent with composition architecture.
    
    properties
        Parameters      % ElegantLaguerreParameters object
        OpticalField    % 2D array representing the field
        r               % radial coordinate matrix
        theta           % angular coordinate matrix
    end
    
    methods
        function obj = ElegantLaguerreBeam(r, theta, params)
            % Constructor
            % r: radial coordinate matrix
            % theta: angular coordinate matrix
            % params: ElegantLaguerreParameters object
            
            if nargin > 0
                obj.Parameters = params;
                obj.r = r;
                obj.theta = theta;
                
                % Elegant argument for Laguerre polynomial: alpha * r^2
                alpha_val = params.alpha;
                xArg = alpha_val .* r.^2;
                
                % Laguerre polynomial: L_p^l(alpha * r^2)
                l = params.l;
                p = params.p;
                Lpl = LaguerreParameters.getAssociatedLaguerrePolynomial(p, l, xArg);
                
                % Physical Gaussian Beam (Carrier)
                % We use the standard GaussianBeam as reference
                GB = GaussianBeam(r, params);
                GField = GB.OpticalField;
                
                % Amplitude term for Elegant beams (standardized):
                % (sqrt(alpha)*r)^|l| 
                % Note: alpha is complex
                amplitudeTerm = (sqrt(alpha_val) .* r).^abs(l);
                
                % Phase and Mode assembly
                % Normalization constant (Simplified for now, as in previous legacy code)
                normConst = 1.0; 
                
                % The final field combines the Gaussian carrier with the elegant polynomial
                obj.OpticalField = normConst .* amplitudeTerm .* Lpl .* exp(1i * l * theta) .* GField;
            end
        end
    end
end
