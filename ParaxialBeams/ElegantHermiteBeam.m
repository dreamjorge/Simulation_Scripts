classdef ElegantHermiteBeam < ParaxialBeam
    % ElegantHermiteBeam - Elegant Hermite-Gaussian beam implementation
    
    properties
        Parameters      % ElegantHermiteParameters object
        OpticalField    % 2D array representing the field
        X               % X coordinate matrix
        Y               % Y coordinate matrix
    end
    
    methods
        function obj = ElegantHermiteBeam(X, Y, params)
            % Constructor
            obj = obj@ParaxialBeam(params.Lambda);
            
            if nargin > 0
                obj.Parameters = params;
                obj.X = X;
                obj.Y = Y;
                obj.OpticalField = obj.computeComplexField(X, Y, params);
            end
        end
        
        function field = opticalField(obj, X, Y, z)
            % Unified API
            field = obj.computeComplexField(X, Y, obj.Parameters);
        end
        
        function field = computeComplexField(obj, x, y, params)
            obj.X = x;
            obj.Y = y;
            
            % Elegant argument for Hermite polynomials: sqrt(alpha) * coordinate
            sqrt_alpha = sqrt(params.alpha);
            argX = sqrt_alpha .* x;
            argY = sqrt_alpha .* y;
            
            % Hermite polynomials: H_n(sqrt(alpha)*x)
            Hn = PolynomialUtils.hermitePoly(params.n, argX);
            Hm = PolynomialUtils.hermitePoly(params.m, argY);
            
            % Physical Gaussian Beam (Carrier)
            [R, ~] = cart2pol(x, y);
            GB = GaussianBeam(R, params);
            GField = GB.OpticalField;
            
            field = Hn .* Hm .* exp(1i * params.PhiPhase) .* GField;
            obj.OpticalField = field;
        end
    end
end
