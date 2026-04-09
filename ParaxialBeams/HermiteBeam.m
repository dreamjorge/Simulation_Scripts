classdef HermiteBeam < ParaxialBeam
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
            obj = obj@ParaxialBeam(params.Lambda);
            
            if nargin > 0
                obj.Parameters = params;
                obj.x = x;
                obj.y = y;
                obj.OpticalField = obj.computeComplexField(x, y, params);
            end
        end
        
        function field = opticalField(obj, X, Y, z)
            % Unified API
            field = obj.computeComplexField(X, Y, obj.Parameters);
        end
        
        function field = computeComplexField(obj, x, y, params)
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
            field = Hn .* Hm .* exp(1i * params.PhiPhase) .* GField;
            obj.OpticalField = field;
        end
    end
    
    methods (Static)
        function H = hermitePoly(n, x)
            H = PolynomialUtils.hermitePoly(n, x);
        end
    end
end