classdef ElegantHermiteBeam
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
            if nargin > 0
                obj.Parameters = params;
                obj.X = X;
                obj.Y = Y;
                
                % Elegant argument for Hermite polynomials: sqrt(alpha) * coordinate
                sqrt_alpha = sqrt(params.alpha);
                argX = sqrt_alpha .* X;
                argY = sqrt_alpha .* Y;
                
                % Hermite polynomials: H_n(sqrt(alpha)*x)
                Hn = HermiteParameters.getHermitePolynomial(params.n, argX);
                Hm = HermiteParameters.getHermitePolynomial(params.m, argY);
                
                % Physical Gaussian Beam (Carrier)
                [R, ~] = cart2pol(X, Y);
                GB = GaussianBeam(R, params);
                GField = GB.OpticalField;
                
                % Normalization constant
                normConst = 1.0;
                
                obj.OpticalField = normConst .* Hn .* Hm .* GField;
            end
        end
    end
end
