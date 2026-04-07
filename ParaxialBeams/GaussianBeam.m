classdef GaussianBeam
    % GaussianBeam - Scalar optical field for a Gaussian beam
    % Compatible with GNU Octave and MATLAB
    
    properties
        Parameters      % GaussianParameters object
        OpticalField    % 2D array representing the field
    end
    
    methods
        function obj = GaussianBeam(r, params)
            % Constructor
            % r: radial coordinate (matrix or array)
            % params: GaussianParameters object
            
            if nargin > 0
                obj.Parameters = params;
                
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
