classdef ElegantLaguerreBeam < ParaxialBeam
    % ElegantLaguerreBeam - Elegant Laguerre-Gaussian beam implementation
    
    properties
        Parameters      % ElegantLaguerreParameters object
        OpticalField    % 2D array representing the field
        r               % radial coordinate matrix
        theta           % angular coordinate matrix
    end
    
    methods
        function obj = ElegantLaguerreBeam(r, theta, params)
            % Constructor
            obj = obj@ParaxialBeam(params.Lambda);
            
            if nargin > 0
                obj.Parameters = params;
                obj.r = r;
                obj.theta = theta;
                obj.OpticalField = obj.computeComplexField(r, theta, params);
            end
        end
        
        function field = opticalField(obj, X, Y, z)
            % Unified API
            [TH, R] = cart2pol(X, Y);
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
            GB = GaussianBeam(r, params);
            GField = GB.OpticalField;
            
            % Amplitude term for Elegant beams
            amplitudeTerm = (sqrt(alpha_val) .* r).^abs(l);
            
            % Mode assembly
            field = amplitudeTerm .* Lpl .* exp(1i * l * theta) .* exp(1i * params.PhiPhase) .* GField;
            obj.OpticalField = field;
        end
    end
end
