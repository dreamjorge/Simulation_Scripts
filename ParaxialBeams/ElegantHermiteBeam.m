classdef ElegantHermiteBeam
    % ElegantHermiteBeam - Elegant Hermite-Gaussian beam implementation
    
    properties
        Parameters      % ElegantHermiteParameters object
        OpticalField    % 2D array representing the field
        x               % x coordinate matrix
        y               % y coordinate matrix
    end
    
    methods
        function obj = ElegantHermiteBeam(x, y, params)
            % Constructor
            if nargin > 0
                obj.Parameters = params;
                obj.x = x;
                obj.y = y;
                
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
                
                obj.OpticalField = Hn .* Hm .* exp(1i * params.PhiPhase) .* GField;
            end
        end
    end
end