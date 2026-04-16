classdef HankelRayTracer < handle
    % HankelRayTracer - Hankel-aware ray propagation with axis-crossing logic.
    %
    % Extends the generic RayTracer approach for Hankel beams (HankelHermite,
    % HankelLaguerre) where each ray carries a Hankel type that may change
    % during propagation.
    %
    % For Laguerre Hankel beams: when a ray crosses the optical axis (r -> 0),
    % H^(2) flips to H^(1) to maintain physical consistency.
    %
    % For Hermite Hankel beams: no axis crossing (Cartesian, no singular axis).
    %
    % Usage:
    %   beam = HankelLaguerre(w0, lambda, l, p, 2);
    %   bundle = RayBundle.createConcentric(Nr, Ntheta, maxR);
    %   bundle.ht(:) = 2;  % all rays start as H^(2)
    %   bundle = HankelRayTracer.propagate(bundle, beam, z_final, dz, 'RK4');

    methods (Static)
        function bundle = propagate(bundle, beam, z_final, dz, method)
            if nargin < 5, method = 'RK4'; end

            z_current = bundle.z(1,1,end);
            while z_current < z_final
                if z_current + dz > z_final
                    dz = z_final - z_current;
                end

                x0 = bundle.x(:,:,end);
                y0 = bundle.y(:,:,end);
                z0 = bundle.z(:,:,end);
                ht0 = bundle.ht(:,:,end);

                r0 = sqrt(x0.^2 + y0.^2);

                if strcmpi(method, 'Euler')
                    [sx, sy] = HankelRayTracer.calculateSlopes(beam, x0, y0, z0, ht0);
                    x1 = x0 + sx .* dz;
                    y1 = y0 + sy .* dz;
                else
                    [k1x, k1y] = HankelRayTracer.calculateSlopes(beam, x0, y0, z0, ht0);
                    [k2x, k2y] = HankelRayTracer.calculateSlopes(beam, x0 + k1x.*dz./2, y0 + k1y.*dz./2, z0 + dz/2, ht0);
                    [k3x, k3y] = HankelRayTracer.calculateSlopes(beam, x0 + k2x.*dz./2, y0 + k2y.*dz./2, z0 + dz/2, ht0);
                    [k4x, k4y] = HankelRayTracer.calculateSlopes(beam, x0 + k3x.*dz, y0 + k3y.*dz, z0 + dz, ht0);

                    sx = (k1x + 2.*k2x + 2.*k3x + k4x) ./ 6;
                    sy = (k1y + 2.*k2y + 2.*k3y + k4y) ./ 6;

                    x1 = x0 + sx .* dz;
                    y1 = y0 + sy .* dz;
                end

                z1 = z0 + dz;

                % Axis-crossing check: if r changed sign, flip H^(2) -> H^(1)
                r1 = sqrt(x1.^2 + y1.^2);
                ht1 = ht0;
                crossed = (r0 > 0) & (r1 < r0) & (r1 < 1e-12);
                ht1(crossed & (ht0 == 2)) = 1;

                bundle.addStep(x1, y1, z1, sx, sy, ht1);
                z_current = z1;
            end
        end

        function [sx, sy] = calculateSlopes(beam, x, y, z, ht)
            delta = 1e-7;
            k = beam.k;

            uniqueTypes = unique(ht(:));
            sx = zeros(size(x));
            sy = zeros(size(x));

            for t = 1:numel(uniqueTypes)
                htype = uniqueTypes(t);
                mask = (ht == htype);

                tempBeam = HankelRayTracer.beamWithType(beam, htype);

                field = tempBeam.opticalField(x, y, z);
                phase = unwrap(angle(field));

                field_dx = tempBeam.opticalField(x + delta, y, z);
                phase_dx = unwrap(angle(field_dx));

                field_dy = tempBeam.opticalField(x, y + delta, z);
                phase_dy = unwrap(angle(field_dy));

                sx_t = (phase_dx - phase) / (delta * k);
                sy_t = (phase_dy - phase) / (delta * k);

                sx(mask) = sx_t(mask);
                sy(mask) = sy_t(mask);
            end
        end
    end

    methods (Static, Access = private)
        function newBeam = beamWithType(beam, htype)
            if isa(beam, 'HankelLaguerre')
                newBeam = HankelLaguerre(beam.InitialWaist, beam.Lambda, beam.l, beam.p, htype);
            elseif isa(beam, 'HankelHermite')
                newBeam = HankelHermite(beam.InitialWaist, beam.Lambda, beam.n, beam.m, htype);
            else
                newBeam = beam;
            end
        end
    end
end
