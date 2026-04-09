classdef GaussianBeam < ParaxialBeam
    % GaussianBeam - Scalar optical field for a Gaussian beam
    % Compatible with GNU Octave and MATLAB
    
    properties
        Parameters      % GaussianParameters object
        OpticalField    % 2D array representing the field (at constructor z)
    end
    
    methods
        function obj = GaussianBeam(r, params)
            % Constructor
            % r: radial coordinate (matrix or array)
            % params: GaussianParameters object
            
            % Call superclass constructor
            obj = obj@ParaxialBeam(params.Lambda);
            
            if nargin > 0
                obj.Parameters = params;
                obj.OpticalField = obj.computeFieldFromR(r, params);
            end
        end
        
        function field = opticalField(obj, X, Y, z)
            % Unified API method
            % X, Y: coordinate matrices
            % z: propagation distance
            
            R_mat = sqrt(X.^2 + Y.^2);
            
            w0 = obj.Parameters.InitialWaist;
            lambda = obj.Lambda;
            zr = pi * w0^2 / lambda;
            
            w = w0 * sqrt(1 + (z ./ zr).^2);
            R = z .* (1 + (zr ./ z).^2);
            % Avoid division by zero at z=0 for R
            if isscalar(z) && z == 0
                R = Inf;
            else
                R(z == 0) = Inf;
            end
            
            psi = atan2(z, zr);
            k = obj.k;
            
            % Amplitude and Gaussian decay
            amplitude = (w0 ./ w) .* exp(-R_mat.^2 ./ w.^2);
            
            % Phase terms
            phase_z = -1i * k .* z;
            phase_curv = 1i * k .* R_mat.^2 ./ (2 .* R);
            phase_curv(isinf(R)) = 0;
            phase_gouy = -1i * psi;
            
            field = amplitude .* exp(phase_z + phase_curv + phase_gouy);
        end
        
        function field = computeFieldFromR(obj, r, params)
             % Internal helper for constructor
             w0 = params.InitialWaist;
             w = params.Waist;
             R = params.Radius;
             k = obj.k;
             z = params.zCoordinate;
             psi = params.GouyPhase;
             
             amplitude = (w0 ./ w) .* exp(-r.^2 ./ w.^2);
             phase_z = -1i * k .* z;
             phase_curv = 1i * k .* r.^2 ./ (2 .* R);
             phase_curv(isinf(R)) = 0;
             phase_gouy = -1i * psi;
             
             field = amplitude .* exp(phase_z + phase_curv + phase_gouy);
        end
    end
end
