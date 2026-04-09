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
            
            % Update/Get parameters for this specific z
            % Note: In a full refactor, we might want to make Parameters dynamic.
            % For now, we use the internal logic.
            R_mat = sqrt(X.^2 + Y.^2);
            
            % We need to calculate parameters at z. 
            % Since GaussianParameters is a value-object/struct-like, 
            % we can instantiate a temporary one or use static methods.
            w0 = obj.Parameters.InitialWaist;
            lambda = obj.Lambda;
            zr = pi * w0^2 / lambda;
            
            w = w0 * sqrt(1 + (z/zr)^2);
            R = z * (1 + (zr/z)^2);
            if z == 0, R = Inf; end
            psi = atan(z/zr);
            k = obj.k;
            
            % Amplitude and Gaussian decay
            amplitude = (w0 / w) .* exp(-R_mat.^2 ./ w.^2);
            
            % Phase terms
            phase_z = -1i * k * z;
            phase_curv = 1i * k * R_mat.^2 ./ (2 * R);
            if isinf(R), phase_curv = zeros(size(R_mat)); end % Handle flat phase at waist
            phase_gouy = -1i * psi;
            
            field = amplitude .* exp(phase_z + phase_curv + phase_gouy);
        end
        
        function field = computeFieldFromR(obj, r, params)
             % Internal helper for constructor
             field = obj.computeField(zeros(size(r)), r, params.zCoordinate);
             % Wait, computeField expects X,Y. If we only have r, we can trick it.
             % Actually, let's just use the logic directly here to avoid confusion.
             w0 = params.InitialWaist;
             w = params.Waist;
             R = params.Radius;
             k = obj.k;
             z = params.zCoordinate;
             psi = params.GouyPhase;
             
             amplitude = (w0 ./ w) .* exp(-r.^2 ./ w.^2);
             phase_z = -1i * k * z;
             phase_curv = 1i * k * r.^2 ./ (2 * R);
             if isinf(R), phase_curv = zeros(size(r)); end
             phase_gouy = -1i * psi;
             
             field = amplitude .* exp(phase_z + phase_curv + phase_gouy);
        end
    end
end
                
                % Calculate complex field
                % Formula: E(r,z) = E0 * (w0/w) * exp(-r^2/w^2) * exp(i*k*z + i*k*r^2/(2*R) - i*psi)
                % Note: Phase convention might vary. Here we follow standard paraxial beam.
                
                w0 = params.InitialWaist;
                w = params.Waist;
                R = params.Radius;
                k = params.k;
                z = params.zCoordinate;
                psi = params.GouyPhase;
                
                % Amplitude and Gaussian decay
                amplitude = (w0 ./ w) .* exp(-r.^2 ./ w.^2);
                
                % Phase terms: propagation, curvature, and Gouy
                % We use -i*k*z convention (consistent with angular spectrum propagator)
                phase_z = -1i * k * z;
                phase_curv = 1i * k * r.^2 ./ (2 * R);
                phase_gouy = -1i * psi;
                
                obj.OpticalField = amplitude .* exp(phase_z + phase_curv + phase_gouy);
            end
        end
    end
end
