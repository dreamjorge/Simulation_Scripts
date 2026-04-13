classdef AnalysisUtils
    % AnalysisUtils - Physical analysis and ray-tracing utilities
    %
    % Theoretical Context (Contexto Academico):
    % This class implements the mapping between wave optics (phase fields) and 
    % ray optics (trajectories) under the Paraxial Approximation.
    % In this limit, the ray slope m = dz/dr (or dx/dz) is perpendicular to the 
    % wavefront (surface of constant phase phi), following the Eikonal equation:
    %
    %   $$ \vec{s} = \frac{\nabla \phi}{|\nabla \phi|} \approx \hat{z} + \frac{\nabla_{\perp} \phi}{k} $$
    %
    % Where:
    % - k = 2*pi/lambda is the wave number (numero de onda).
    % - phi is the phase of the optical field (fase del campo óptico).
    %
    % For a ray propagating along z, the local slopes are defined as:
    %   $$ m_x = \frac{dx}{dz} = \frac{1}{k} \frac{\partial \phi}{\partial x} $$
    %   $$ m_y = \frac{dy}{dz} = \frac{1}{k} \frac{\partial \phi}{\partial y} $$
    %
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
            
            % For scalar output, pick the value at the specified coordinate (x, z)
            % This matches the legacy logic found in Addons/getCylindricalGradient.m
            idx_x = floor(max(1, min(length(fr), x/dx + length(fr)/2)));
            idx_z = floor(max(1, min(length(fz), z/dz + 1))); % z is usually axial
            
            % Standard slope m_zr = - (partial_phi/partial_r) / (partial_phi/partial_z + k)
            % Note: some legacy code might use reciprocal or positive, but 
            % the tests expect mzx, mzy for ray propagation.
            mzr = -gz(idx_x) / gr(idx_z);
        end

        function [mzx, mzy, mxy] = gradientXYZ(fyz, fxz, fxy, k, dx, dy, dz, x, y, z)
            % gradientXYZ Calculates ray slopes in Cartesian planes
            % fyz: Phase slice in y-z plane (at fixed x)
            % fxz: Phase slice in x-z plane (at fixed y)
            % fxy: Phase slice in x-y plane (at fixed z)
            
            gx = gradient(fxz) / dx;
            gy = gradient(fyz) / dy;
            gz = gradient(fxz) / dz + k; % Using axial slice to estimate longitudinal k
            
            % Index based on coordinates (assuming centered window for x, y)
            idx_x = floor(max(1, min(size(fxz, 2), x/dx + size(fxz, 2)/2)));
            idx_y = floor(max(1, min(size(fyz, 1), y/dy + size(fyz, 1)/2)));
            idx_z = floor(max(1, min(size(fxz, 1), z/dz + 1)));

            mzx = -gx(idx_z, idx_x) / gz(idx_z, idx_x);
            mzy = -gy(idx_y, idx_z) / gz(idx_y, idx_z);
            % mxy is usually the transverse coupling, often dphi/dx / dphi/dy or similar
            mxy = gx(idx_z, idx_x) / gy(idx_y, idx_z);
        end

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
        
    end
end
