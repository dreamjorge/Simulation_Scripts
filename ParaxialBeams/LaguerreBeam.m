classdef LaguerreBeam
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
            % r: radial coordinate matrix
            % theta: angular coordinate matrix
            % params: LaguerreParameters object
            
            if nargin > 0
                obj.Parameters = params;
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

                % GField = (w0/w)*exp(-r^2/w^2)*exp(i*phase). Gouy for LG adds (2p+|l|)*psi via PhiPhase.
                obj.OpticalField = amplitudeTerm .* Lpl .* exp(1i * l * theta) .* exp(1i * params.PhiPhase) .* GField;
            end
        end
    end
end