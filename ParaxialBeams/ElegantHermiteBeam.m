classdef ElegantHermiteBeam
    % ElegantHermiteBeam - Elegant Hermite-Gaussian beam implementation
    %
    % Coordinate system: CARTESIAN (x, y) - same as HermiteBeam.
    %
    % Formula differences from standard HermiteBeam:
    %   Standard HG: H_n(sqrt(2)*x/w), H_m(sqrt(2)*y/w)
    %   Elegant HG: H_n(sqrt(alpha)*x), H_m(sqrt(alpha)*y)
    %
    % The "elegant" variant uses the complex beam parameter alpha (not the waist w)
    % to define the polynomial scaling. This is the "elegant" Hermite-Gauss convention
    % from Siegman (1990).
    
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
                R = sqrt(x.^2 + y.^2);
                GB = GaussianBeam(R, params);
                GField = GB.OpticalField;
                
                obj.OpticalField = Hn .* Hm .* exp(1i * params.PhiPhase) .* GField;
            end
        end
    end
end