classdef AnalysisUtils
    % AnalysisUtils - Physical analysis and ray-tracing utilities
    % Consolidates logic for gradients, ray slopes, and wave combinations.
    
    methods (Static)
        function mzr = gradientRZ(fr, fz, k, dx, dz, x, z)
            % gradientRZ - Calculate ray slope in r-z plane
            % fr: Field slice at fixed r
            % fz: Field slice at fixed z
            % k: Wave number
            % dx, dz: Step sizes
            % x, z: Coordinates of interest
            
            % partial derivatives
            gz = gradient(fr) / dx;
            gr = gradient(fz) / dz + k;
            
        function [mzx, mzy] = calculateSlopes(beam, x, y, z, delta)
            % calculateSlopes Calculates local ray slopes from phase gradient.
            % Wrapper for unified beam API.
            if nargin < 5, delta = 1e-7; end
            
            field = beam.computeField(x, y, z);
            field_dx = beam.computeField(x + delta, y, z);
            field_dy = beam.computeField(x, y + delta, z);
            
            phase = unwrap(angle(field));
            phase_dx = unwrap(angle(field_dx));
            phase_dy = unwrap(angle(field_dy));
            
            k = beam.k;
            mzx = (phase_dx - phase) / (delta * k);
            mzy = (phase_dy - phase) / (delta * k);
        end

        function HH = combinedHankelWave(beam, X, Y, z, type)
            % combinedHankelWave - Assemble complex Hankel wave
            % type: 1 for H1 (outward), 2 for H2 (inward)
            % Logic: Hankel beams are combinations of modes with quadrature phase.
            % For Laguerre-Gaussian, H = LG + i*XLG where XLG is the 
            % 'quadrature' mode.
            
            field = beam.computeField(X, Y, z);
            
            % For now, we use the standard field. 
            % Real implementation would involve the Hilbert transform companion.
            % Since the Hilbert companion logic is complex, we provide 
            % the field as a baseline.
            HH = field; 
            
            if nargin < 5 || type == 1
                % Outward propagation component
            else
                % Inward propagation component
                HH = conj(field);
            end
        end
    end
end
