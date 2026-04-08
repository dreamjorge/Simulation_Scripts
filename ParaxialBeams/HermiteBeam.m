classdef HermiteBeam
    % HermiteBeam - Scalar optical field for a Hermite-Gaussian beam
    % Compatible with GNU Octave and MATLAB
    
    properties
        Parameters      % HermiteParameters object
        OpticalField    % 2D array representing the field
        x               % x coordinate matrix
        y               % y coordinate matrix
    end
    
    methods
        function obj = HermiteBeam(x, y, params)
            % Constructor
            % x, y: coordinate matrices
            % params: HermiteParameters object
            
            if nargin > 0
                obj.Parameters = params;
                obj.x = x;
                obj.y = y;
                
                % Radial coordinate for Gaussian part
                r = sqrt(x.^2 + y.^2);
                
                % Fundametal Gaussian Field
                GB = GaussianBeam(r, params);
                GField = GB.OpticalField;
                
                % Hermite polynomials part
                w = params.Waist;
                n = params.n;
                m = params.m;
                
                Hn = PolynomialUtils.hermitePoly(n, sqrt(2) * x ./ w);
                Hm = PolynomialUtils.hermitePoly(m, sqrt(2) * y ./ w);
                
                % Phase shift (n+m)*psi is handled in Parameters.PhiPhase
                % We multiply by the phase term and the polynomials
                obj.OpticalField = Hn .* Hm .* exp(1i * params.PhiPhase) .* GField;
            end
        end
    end
    
    methods (Static)
        function H = hermitePoly(n, x)
            % Kept for backwards compatibility. Delegates to PolynomialUtils.
            H = PolynomialUtils.hermitePoly(n, x);
        end
    end
end