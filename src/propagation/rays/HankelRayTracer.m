classdef HankelRayTracer < handle
    % HankelRayTracer — Phase-gradient ray tracing for Hankel (cylindrical) beams.
    %
    % ============================================================================
    % HANKEL BEAMS AND THEIR PHYSICS
    % ============================================================================
    %
    % Hankel beams are solutions of the paraxial wave equation in cylindrical
    % coordinates (r, θ, z). Unlike Hermite-Gaussian beams (Cartesian), they
    % carry orbital angular momentum (OAM) and have a phase singularity on axis.
    %
    % The complex amplitude of a Hankel-Laguerre beam is:
    %
    %   u(r, θ, z) = w(z)/w₀ · (√2·r/w(z))^{|l|}
    %                · exp(-r²/w(z)²) · exp(i·l·θ)
    %                · exp(-i·k·z - i·k·r²/(2·R(z)) + i·ψ(z))
    %
    % where:
    %   l    = azimuthal mode index (integer, OAM per photon = l·ℏ)
    %   p    = radial mode index (integer)
    %   w(z) = beam radius at z
    %   R(z) = wavefront radius of curvature
    %   ψ(z) = Gouy phase = arctan(z/zᵣ)
    %
    % At r = 0, the amplitude vanishes for l ≠ 0 (phase singularity/vortex).
    % The phase is undefined at the axis, and the field has a screw dislocation.
    %
    % Hankel beams are described in:
    %   - Allen92: L. Allen et al., "Orbital angular momentum of
    %     light and the transformation of Laguerre-Gaussian laser modes",
    %     Phys. Rev. A 45, 1992.
    %   - Papi: J. Papi, "Hankel beams", in "Structured Light", Elsevier 2023.
    %
    % ============================================================================
    % BRANCH SWITCHING: WHY H^(1) AND H^(2) MATTER
    % ============================================================================
    %
    % The radial part of Hankel beams is expressed as a superposition of
    % inward (H^(1)) and outward (H^(2)) cylindrical waves:
    %
    %   H_l^(1)(r,z) ~ J_l(kr) + i·Y_l(kr)   (inward propagating)
    %   H_l^(2)(r,z) ~ J_l(kr) - i·Y_l(kr)   (outward propagating)
    %
    % At the axis (r → 0), both J_l and Y_l behave as:
    %   J_l(r) ~ (r/2)^{|l|} / Γ(|l|+1)   (finite for l≠0)
    %   Y_l(r) ~ -Γ(|l|)·(2/r)^{|l|}      (singular for l≠0)
    %
    % The combination H_l^(1) or H_l^(2) determines the vorticity direction
    % (phase winding sense around axis). When a ray crosses the axis, the
    % Hankel branch must flip to maintain physical continuity — the inward
    % wave becomes outward and vice versa.
    %
    % The flip condition is geometric: when the ray segment passes within
    % a small distance of the axis, we switch the branch to maintain the
    % correct topological charge.
    %
    % References:
    %   - Nienhuis96: G. Nienhuis, "Doppler effect induced by rotating
    %     lenses", Opt. Commun. 119 (1996).
    %   - Courtial98: J. Courtial et al., "Gaussian beams with
    %     singular phase", Pure Appl. Opt. 7 (1998).
    %
    % ============================================================================
    % AXIS-CROSSING DETECTION: MINIMUM SEGMENT DISTANCE
    % ============================================================================
    %
    % Previous implementation used the sign of the 2D cross product:
    %   det = x₀·y₁ - x₁·y₀ = |r₀|·|r₁|·sin(θ)
    % as a proxy for "did the ray cross the axis". This is WRONG:
    %
    %   - det only measures the orientation (signed area of triangle O,r₀,r₁).
    %   - det < 0 can occur when the ray curves near the axis without
    %     actually passing through it (grazing trajectory).
    %   - det changes sign when the ray goes "behind" the axis from the
    %     viewpoint, even if it doesn't physically cross.
    %
    % Correct criterion: minimum Euclidean distance from the ray segment
    % to the origin.
    %
    % For segment P(s) = P₀ + s·(P₁-P₀) with s ∈ [0,1]:
    %   The closest point to origin occurs at s* = -P₀·(P₁-P₀) / |P₁-P₀|²
    %   The minimum distance is |P₀ + s*·(P₁-P₀)|.
    %
    % We flip when min_dist < threshold, where the threshold adapts to
    % the local geometry: threshold = max(|x₀|, |y₀|) · 10⁻³.
    %
    % This ensures:
    %   (a) Rays aimed at origin (small impact parameter) flip.
    %   (b) Grazing rays that don't actually cross don't flip.
    %   (c) The flip point is geometrically well-defined.
    %
    % References:
    %   - Wikipedia: "Distance from a point to a line" (segment case).
    %   - Eric & Kucnera: "Minimum distance calculation for line segments",
    %     J. Graphics Tools 8, 2003.
    %
    % ============================================================================

    methods (Static)

        function bundle = propagate(bundle, beam, z_final, dz, method)
            % PROPAGATE — Integrate Hankel ray bundle with axis-crossing detection.
            %
            % Same as RayTracer.propagate but tracks the Hankel type (ht)
            % per ray. When a ray segment passes within threshold distance
            % of the optical axis, ht flips from 2 → 1 (outward ↔ inward).
            %
            % See class documentation for axis-crossing criterion.

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

                if strcmpi(method, 'Euler')
                    [sx, sy] = HankelRayTracer.calculateSlopes(beam, x0, y0, z0, ht0);
                    x1 = x0 + sx .* dz;
                    y1 = y0 + sy .* dz;
                else
                    [k1x, k1y] = HankelRayTracer.calculateSlopes(beam, x0,           y0,           z0, ht0);
                    [k2x, k2y] = HankelRayTracer.calculateSlopes(beam, x0+k1x.*dz/2, y0+k1y.*dz/2, z0+dz/2, ht0);
                    [k3x, k3y] = HankelRayTracer.calculateSlopes(beam, x0+k2x.*dz/2, y0+k2y.*dz/2, z0+dz/2, ht0);
                    [k4x, k4y] = HankelRayTracer.calculateSlopes(beam, x0+k3x.*dz,   y0+k3y.*dz,   z0+dz,   ht0);

                    sx = (k1x + 2.*k2x + 2.*k3x + k4x) ./ 6;
                    sy = (k1y + 2.*k2y + 2.*k3y + k4y) ./ 6;

                    x1 = x0 + sx .* dz;
                    y1 = y0 + sy .* dz;
                end

                z1 = z0 + dz;

                % ----------------------------------------------------------------
                % AXIS-CROSSING CHECK
                %
                % Given segment from P₀ = (x₀,y₀) to P₁ = (x₁,y₁),
                % find the point on the segment closest to origin.
                %
                % Parametric form: P(s) = P₀ + s·(P₁-P₀), s ∈ [0,1]
                %
                % Distance to origin squared: |P(s)|²
                % Derivative w.r.t. s: d/ds|P(s)|² = 2·P₀·(P₁-P₀) + 2·s·|P₁-P₀|²
                %
                % Setting to zero:
                %   s* = -P₀·(P₁-P₀) / |P₁-P₀|²
                %
                % Clamp s* to [0,1] → projection falls within segment.
                % Minimum distance = |P₀ + s*_clamped·(P₁-P₀)|
                %
                % Threshold = max(|x₀|,|y₀|) · 10⁻³ adapts to local scale.
                % ----------------------------------------------------------------
                dx = x1 - x0;
                dy = y1 - y0;
                denominator = dx.^2 + dy.^2 + eps;

                % s* = -P₀·(P₁-P₀) / |P₁-P₀|²
                t_param = -(x0 .* dx + y0 .* dy) ./ denominator;

                % Clamp projection to segment bounds
                t_clamped = max(0, min(1, t_param));

                % Nearest point on segment to origin
                nearest_x = x0 + t_clamped .* dx;
                nearest_y = y0 + t_clamped .* dy;

                % Minimum Euclidean distance
                min_dist = sqrt(nearest_x.^2 + nearest_y.^2);

                % Adaptive threshold: local scale × tolerance
                threshold = max(abs(x0), abs(y0)) * 1e-3;

                % Flip only if segment passes close enough to axis
                crossed = min_dist < threshold;

                ht1 = ht0;
                ht1(crossed & (ht0 == 2)) = 1;  % flip H^(2) → H^(1)

                bundle.addStep(x1, y1, z1, sx, sy, ht1);
                z_current = z1;
            end
        end


        function [sx, sy] = calculateSlopes(beam, x, y, z, ht)
            % CALCULATESLOPES — Phase gradient for Hankel beams by type.
            %
            % Hankel beams require different field evaluations for each
            % branch (H^(1) vs H^(2)). We group rays by their current
            % ht value and evaluate the corresponding field for each group.
            %
            % The gradient method is identical to RayTracer
            % (complex field phase gradient), but applied per-type.

            epsilon = 1e-12;
            w0     = beam.InitialWaist;
            lambda = beam.Lambda;
            delta  = RayTracer.resolveDelta(x, y, w0, lambda);

            uniqueTypes = unique(ht(:));
            sx = zeros(size(x));
            sy = zeros(size(x));

            for t = 1:numel(uniqueTypes)
                htype = uniqueTypes(t);
                mask  = (ht == htype);

                % Build beam object with correct Hankel branch
                tempBeam = HankelRayTracer.beamWithType(beam, htype);

                % Complex field phase gradient (same as RayTracer)
                u0   = tempBeam.opticalField(x, y, z);
                u_xp = tempBeam.opticalField(x+delta, y, z);
                u_xm = tempBeam.opticalField(x-delta, y, z);
                u_yp = tempBeam.opticalField(x, y+delta, z);
                u_ym = tempBeam.opticalField(x, y-delta, z);

                dudx = (u_xp - u_xm) / (2 * delta);
                dudy = (u_yp - u_ym) / (2 * delta);

                u0_conj  = conj(u0);
                abs_u0_sq = real(u0_conj .* u0);

                sx_num = imag(u0_conj .* dudx);
                sy_num = imag(u0_conj .* dudy);

                denom = abs_u0_sq + epsilon;
                sx_t = sx_num ./ denom;
                sy_t = sy_num ./ denom;

                sx(mask) = sx_t(mask);
                sy(mask) = sy_t(mask);
            end
        end

    end % methods (Static)


    methods (Static, Access = private)

        function newBeam = beamWithType(beam, htype)
            % BEAMWITHTYPE — Create a beam copy with a specific Hankel branch.
            %
            % Hankel beams carry a type index: 1 = H^(1) (inward),
            % 2 = H^(2) (outward). This method fabricates the
            % correct beam variant without re-computing the field formula.
            %
            % For HankelLaguerre: the Hankel type determines whether
            % the radial Bessel function uses H^(1) or H^(2), which
            % have different asymptotic behavior at infinity.
            %
            % For Hermite beams: no branch switching (Cartesian,
            % no singular axis) — returns original beam unchanged.

            if isa(beam, 'HankelLaguerre')
                newBeam = HankelLaguerre(beam.InitialWaist, beam.Lambda, ...
                    beam.l, beam.p, htype);
            elseif isa(beam, 'HankelHermite')
                newBeam = HankelHermite(beam.InitialWaist, beam.Lambda, ...
                    beam.n, beam.m, htype);
            else
                newBeam = beam;
            end
        end

    end % methods (Static, Access = private)

end % classdef