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
            % Uses complex field method: ∇φ = Im(u̅∇u) / |u|²
            % Falls back to central-difference phase gradient if direct method fails.
            sx = 0; sy = 0; % ensure output vars exist on early return
            try
                [sx, sy] = RayTracer.calculatePhaseGradientComplex(beam, x, y, z);
            catch
                % Fallback to central-difference phase gradient
                w0 = beam.InitialWaist;
                lambda = beam.Lambda;
                delta = RayTracer.resolveDelta(x, y, w0, lambda);
                k = beam.k;
                field = beam.opticalField(x, y, z);
                phase = unwrap(angle(field));
                field_dx = beam.opticalField(x + delta, y, z);
                phase_dx = unwrap(angle(field_dx));
                field_dy = beam.opticalField(x, y + delta, z);
                phase_dy = unwrap(angle(field_dy));
                sx = (phase_dx - phase) / (delta * k);
                sy = (phase_dy - phase) / (delta * k);
            end
        end

        function [sx, sy] = calculatePhaseGradientComplex(beam, x, y, z)
            % Compute phase gradient: ∇φ = Im(u̅∇u) / (|u|² + ε)
            % Uses central difference for spatial derivatives of field.
            % Returns sx, sy with same units as dx/dz, dy/dz (radians/m).
            epsilon = 1e-12;
            w0 = beam.InitialWaist;
            lambda = beam.Lambda;
            delta = RayTracer.resolveDelta(x, y, w0, lambda);

            % Field at center and offsets for central difference
            u0  = beam.opticalField(x, y, z);
            u_xp = beam.opticalField(x + delta, y, z);
            u_xm = beam.opticalField(x - delta, y, z);
            u_yp = beam.opticalField(x, y + delta, z);
            u_ym = beam.opticalField(x, y - delta, z);

            % Central difference for ∂u/∂x and ∂u/∂y
            dudx = (u_xp - u_xm) / (2 * delta);
            dudy = (u_yp - u_ym) / (2 * delta);

            % Phase gradient via complex field identity:
            % ∇φ = Im(conj(u) * ∇u) / |u|²
            u0_conj = conj(u0);
            abs_u0_sq = real(u0_conj .* u0);  % |u|²

            sx_num = imag(u0_conj .* dudx);
            sy_num = imag(u0_conj .* dudy);

            denominator = abs_u0_sq + epsilon;
            sx = sx_num ./ denominator;
            sy = sy_num ./ denominator;
        end

        function delta = resolveDelta(x, y, w0, lambda)
            % Adaptive delta based on wavelength and local spatial scale.
            % Prevents cancellation at small x,y while staying large enough
            % to capture field curvature at large radii.
            delta = max(lambda, abs(x) * 1e-4, abs(y) * 1e-4, w0 * 1e-4);
        end
    end
end
