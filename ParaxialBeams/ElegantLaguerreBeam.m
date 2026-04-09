classdef ElegantLaguerreBeam < ParaxialBeam
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
            
            obj = obj@ParaxialBeam(params.Lambda);
            
            if nargin > 0
                obj.Parameters = params;
                obj.r = r;
                obj.theta = theta;
                obj.OpticalField = obj.computeComplexField(r, theta, params);
            end
        end
        
        function field = opticalField(obj, X, Y, z)
            [TH, R] = cart2pol(X, Y);
            % Unified API
            field = obj.computeComplexField(R, TH, obj.Parameters);
        end
        
        function field = computeComplexField(obj, r, theta, params)
                obj.r = r;
                obj.theta = theta;
                
                % Elegant argument for Laguerre polynomial: alpha * r^2
                alpha_val = params.alpha;
                xArg = alpha_val .* r.^2;
                
                % Laguerre polynomial: L_p^l(alpha * r^2)
                l = params.l;
                p = params.p;
                Lpl = PolynomialUtils.associatedLaguerre(p, l, xArg);
                
                % Physical Gaussian Beam (Carrier)
                % We use the standard GaussianBeam as reference
                GB = GaussianBeam(r, params);
                GField = GB.OpticalField;
                
                % Amplitude term for Elegant beams (standardized):
                % (sqrt(alpha)*r)^|l| 
                % Note: alpha is complex
                amplitudeTerm = (sqrt(alpha_val) .* r).^abs(l);
                
                % Phase and mode assembly (Gouy shift consistent with standard LG beam)
                obj.OpticalField = amplitudeTerm .* Lpl .* exp(1i * l * theta) .* exp(1i * params.PhiPhase) .* GField;
            end
        end
    end
end
