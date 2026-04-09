classdef RayTracer < handle
    % RayTracer Utility class for optical ray propagation.
    
    methods (Static)
        function bundle = propagate(bundle, beam, z_final, dz, method)
            % bundle: RayBundle object
            % beam: Beam object (GaussianBeam, etc.)
            % z_final: Final z coordinate
            % dz: Step size
            % method: 'Euler' or 'RK4'
            
            if nargin < 5, method = 'RK4'; end
            
            z_current = bundle.z(1,1,end);
            while z_current < z_final
                % Ensure we don't overshoot
                if z_current + dz > z_final
                    dz = z_final - z_current;
                end
                
                % Get current state (last slice)
                x0 = bundle.x(:,:,end);
                y0 = bundle.y(:,:,end);
                z0 = bundle.z(:,:,end);
                
                if strcmpi(method, 'Euler')
                    [sx, sy] = RayTracer.calculateSlopes(beam, x0, y0, z0);
                    x1 = x0 + sx * dz;
                    y1 = y0 + sy * dz;
                else
                    % RK4
                    [k1x, k1y] = RayTracer.calculateSlopes(beam, x0, y0, z0);
                    [k2x, k2y] = RayTracer.calculateSlopes(beam, x0 + k1x*dz/2, y0 + k1y*dz/2, z0 + dz/2);
                    [k3x, k3y] = RayTracer.calculateSlopes(beam, x0 + k2x*dz/2, y0 + k2y*dz/2, z0 + dz/2);
                    [k4x, k4y] = RayTracer.calculateSlopes(beam, x0 + k3x*dz, y0 + k3y*dz, z0 + dz);
                    
                    sx = (k1x + 2*k2x + 2*k3x + k4x) / 6;
                    sy = (k1y + 2*k2y + 2*k3y + k4y) / 6;
                    
                    x1 = x0 + sx * dz;
                    y1 = y0 + sy * dz;
                end
                
                z1 = z0 + dz;
                bundle.addStep(x1, y1, z1, sx, sy);
                z_current = z1;
            end
        end
        
        function [sx, sy] = calculateSlopes(beam, x, y, z)
            % Calculate local phase gradients as slopes dx/dz, dy/dz
            % We use the analytical phase if available, otherwise numerical
            
            % For paraxial beams: dx/dz = (1/k) * d(phase)/dx
            % In our implementation, we'll use a small finite difference for gradient
            % because not all beams might have analytical derivatives implemented.
            
            delta = 1e-7; % Small perturbation for numerical gradient
            
            field = beam.opticalField(x, y, z);
            phase = unwrap(angle(field));
            
            field_dx = beam.opticalField(x + delta, y, z);
            phase_dx = unwrap(angle(field_dx));
            
            field_dy = beam.opticalField(x, y + delta, z);
            phase_dy = unwrap(angle(field_dy));
            
            k = 2*pi / beam.Lambda;
            
            sx = (phase_dx - phase) / (delta * k);
            sy = (phase_dy - phase) / (delta * k);
        end
    end
end
