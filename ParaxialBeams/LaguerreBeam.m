classdef LaguerreBeam < ParaxialBeam
    % LaguerreBeam - Scalar optical field for a Laguerre-Gaussian beam
    % Compatible with GNU Octave and MATLAB
    %
    % Coordinate system: POLAR (r, theta)
    % Note: This differs from HermiteBeam which uses Cartesian (x, y).
    % The beam classes do not share a unified coordinate API.
    
    properties
        Parameters      % LaguerreParameters object
        OpticalField    % 2D array representing the field
        r               % radial coordinate matrix
        theta           % angular coordinate matrix
    end
    
    methods
        function obj = LaguerreBeam(r, theta, params)
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
            
            % Fundamental Gaussian Field
            GB = GaussianBeam(r, params);
            GField = GB.OpticalField;
            
            % Laguerre mode part
            w = params.Waist;
            l = params.l;
            p = params.p;
            
            % Amplitude term: (sqrt(2)*r/w)^|l|
            amplitudeTerm = (sqrt(2) * r ./ w).^abs(l);

            % Laguerre polynomial: L_p^|l|(2*r^2 / w^2)
            xArg = 2 * r.^2 ./ w.^2;
            Lpl = PolynomialUtils.associatedLaguerre(p, l, xArg);

            % Mode assembly
            field = amplitudeTerm .* Lpl .* exp(1i * l * theta) .* exp(1i * params.PhiPhase) .* GField;
            obj.OpticalField = field;
        end
    end
end