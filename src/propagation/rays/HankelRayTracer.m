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

        function bundleOut = propagateToPlanes(bundleIn, beam, zPlanes, dzInternal, method)
            % PROPAGATETOPLANES — Propagate with internal substeps, sample at fixed z planes.
            %
            % Why this exists:
            %   Field propagation/visualization is usually rendered on fixed z planes
            %   (e.g. FFT planes). If ray integration changes internal step size,
            %   indexing by step can desynchronize rays vs images.
            %
            % This method keeps overlay-consistent outputs by always returning ray
            % states exactly at user-specified z planes, while using dzInternal for
            % internal integration between planes.
            %
            % Inputs:
            %   bundleIn   : initial RayBundle (uses last slice as starting state)
            %   beam       : beam model
            %   zPlanes    : monotonically increasing vector of z sample planes
            %   dzInternal : internal integration step (scalar > 0)
            %   method     : 'RK4' (default) or 'Euler'
            %
            % Output:
            %   bundleOut  : RayBundle sampled exactly at zPlanes

            if nargin < 5 || isempty(method), method = 'RK4'; end
            if nargin < 4 || isempty(dzInternal), dzInternal = []; end

            if isempty(zPlanes)
                bundleOut = bundleIn;
                return;
            end

            zPlanes = zPlanes(:).';
            if any(diff(zPlanes) < 0)
                error('zPlanes must be monotonically increasing.');
            end

            % Current state from input bundle last slice
            x0 = bundleIn.x(:,:,end);
            y0 = bundleIn.y(:,:,end);
            z0 = bundleIn.z(:,:,end);
            ht0 = bundleIn.ht(:,:,end);

            zStart = z0(1,1);
            if abs(zPlanes(1) - zStart) > 1e-15
                error('zPlanes(1) must match initial bundle z (%.6e).', zStart);
            end

            % Output bundle starts at initial state only (fixed-plane samples)
            bundleOut = RayBundle(x0, y0, zStart);
            bundleOut.ht = ht0;

            for kk = 2:numel(zPlanes)
                zTarget = zPlanes(kk);
                if zTarget < zStart
                    error('zPlanes must be increasing and >= initial z.');
                end

                % Internal step for this interval
                if isempty(dzInternal)
                    dzStep = max((zTarget - zStart) / 20, eps);
                else
                    dzStep = dzInternal;
                end

                % Propagate from current state to current target using internal steps
                bTmp = RayBundle(x0, y0, zStart);
                bTmp.ht = ht0;
                bTmp = HankelRayTracer.propagate(bTmp, beam, zTarget, dzStep, method);

                x1 = bTmp.x(:,:,end);
                y1 = bTmp.y(:,:,end);
                z1 = bTmp.z(:,:,end);
                sx1 = bTmp.sx(:,:,end);
                sy1 = bTmp.sy(:,:,end);
                ht1 = bTmp.ht(:,:,end);

                bundleOut.addStep(x1, y1, z1, sx1, sy1, ht1);

                % advance state
                x0 = x1; y0 = y1; z0 = z1; ht0 = ht1;
                zStart = z0(1,1);
            end
        end

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

                % ----------------------------------------------------------------
                % SUB-STEP CROSSING CORRECTION (smooth H2->H1 transition)
                %
                % Problem: flipping only at end-of-step causes a visible kink near
                % origin (piecewise slope jump) when a ray crosses the axis.
                %
                % Fix: for crossed rays, split integration into two sub-steps:
                %   1) integrate with current branch up to crossing fraction t*
                %   2) switch branch (2->1) and integrate remaining segment
                %
                % This keeps position update smooth through the crossing point and
                % matches legacy behavior where H^(2) rays pass origin then invert.
                % ----------------------------------------------------------------
                if any(crossed(:))
                    idxCross = find(crossed);
                    for ii = 1:numel(idxCross)
                        idx = idxCross(ii);

                        x0i = x0(idx); y0i = y0(idx); z0i = z0(idx);
                        ht0i = ht0(idx);

                        % Crossing fraction along this step (already clamped [0,1])
                        tCross = t_clamped(idx);
                        dz1 = dz * tCross;
                        dz2 = dz - dz1;

                        % Step 1: up to crossing with original branch
                        xMid = x0i; yMid = y0i; zMid = z0i;
                        if dz1 > 0
                            [xMid, yMid] = HankelRayTracer.singleStep(beam, x0i, y0i, z0i, ht0i, dz1, method);
                            zMid = z0i + dz1;
                        end

                        % Step 2: after crossing, force H^(2)->H^(1)
                        htMid = ht0i;
                        if ht0i == 2
                            htMid = 1;
                        end

                        x1i = xMid; y1i = yMid;
                        if dz2 > 0
                            [x1i, y1i] = HankelRayTracer.singleStep(beam, xMid, yMid, zMid, htMid, dz2, method);
                        end

                        x1(idx) = x1i;
                        y1(idx) = y1i;
                        ht1(idx) = htMid;

                        % Effective slope stored for this global step
                        sx(idx) = (x1i - x0i) / dz;
                        sy(idx) = (y1i - y0i) / dz;
                    end
                end

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
            % For beams WITH vortex (HankelLaguerre with l ≠ 0), we use
            % polar-coordinate gradient (calculatePhaseGradientPolar) to avoid
            % spurious gradients from branch-cut crossing at θ = 0, π.
            %
            % For beams WITHOUT vortex, we use the same Cartesian complex
            % gradient as RayTracer.

            uniqueTypes = unique(ht(:));
            sx = zeros(size(x));
            sy = zeros(size(x));

            % ----------------------------------------------------------------
            % Smooth H^(2)->H^(1) transition for non-vortex Laguerre (l=0)
            %
            % For l=0 there is no azimuthal singularity. A hard branch switch
            % close to axis can still introduce a visible kink in trajectories.
            % Blend slopes smoothly near origin for rays currently in H^(2).
            % ----------------------------------------------------------------
            if isa(beam, 'HankelLaguerre') && (beam.l == 0)
                b1 = HankelRayTracer.beamWithType(beam, 1);
                b2 = HankelRayTracer.beamWithType(beam, 2);

                [sx1, sy1] = RayTracer.calculatePhaseGradientComplex(b1, x, y, z);
                [sx2, sy2] = RayTracer.calculatePhaseGradientComplex(b2, x, y, z);

                % Transition radius: small core around axis
                r = sqrt(x.^2 + y.^2);
                rBlend = max(8 * beam.Lambda, 0.06 * beam.InitialWaist);
                rWidth = max(2 * beam.Lambda, 0.25 * rBlend);
                alpha = 0.5 * (1 - tanh((r - rBlend) ./ max(rWidth, eps)));

                % H^(1) rays keep H1 slope; H^(2) rays smoothly morph to H1 near axis
                mask1 = (ht == 1);
                mask2 = (ht == 2);

                sxMix = (1 - alpha) .* sx2 + alpha .* sx1;
                syMix = (1 - alpha) .* sy2 + alpha .* sy1;

                sx(mask1) = sx1(mask1);
                sy(mask1) = sy1(mask1);
                sx(mask2) = sxMix(mask2);
                sy(mask2) = syMix(mask2);
                return;
            end

            if RayTracer.beamHasVortex(beam)
                % Vortex beam: use polar-coordinate gradient
                for t = 1:numel(uniqueTypes)
                    htype = uniqueTypes(t);
                    mask  = (ht == htype);

                    % Build beam object with correct Hankel branch
                    tempBeam = HankelRayTracer.beamWithType(beam, htype);

                    % Polar gradient (handles vortex naturally)
                    [sx_part, sy_part] = RayTracer.calculatePhaseGradientPolar(tempBeam, x, y, z);

                    sx(mask) = sx_part(mask);
                    sy(mask) = sy_part(mask);
                end
            else
                for t = 1:numel(uniqueTypes)
                    htype = uniqueTypes(t);
                    mask  = (ht == htype);

                    tempBeam = HankelRayTracer.beamWithType(beam, htype);
                    [sx_part, sy_part] = RayTracer.calculatePhaseGradientComplex(tempBeam, x, y, z);

                    sx(mask) = sx_part(mask);
                    sy(mask) = sy_part(mask);
                end
            end
        end

    end % methods (Static)


    methods (Static, Access = private)

        function [x1, y1] = singleStep(beam, x0, y0, z0, ht0, dz, method)
            % SINGLESTEP — Integrate one ray over a scalar dz.
            % Used by crossing sub-step correction to split H2->H1 transition.

            if strcmpi(method, 'Euler')
                [sx, sy] = HankelRayTracer.calculateSlopes(beam, x0, y0, z0, ht0);
                x1 = x0 + sx .* dz;
                y1 = y0 + sy .* dz;
            else
                [k1x, k1y] = HankelRayTracer.calculateSlopes(beam, x0,              y0,              z0,       ht0);
                [k2x, k2y] = HankelRayTracer.calculateSlopes(beam, x0+k1x.*dz/2,    y0+k1y.*dz/2,    z0+dz/2,  ht0);
                [k3x, k3y] = HankelRayTracer.calculateSlopes(beam, x0+k2x.*dz/2,    y0+k2y.*dz/2,    z0+dz/2,  ht0);
                [k4x, k4y] = HankelRayTracer.calculateSlopes(beam, x0+k3x.*dz,      y0+k3y.*dz,      z0+dz,    ht0);

                sx = (k1x + 2.*k2x + 2.*k3x + k4x) ./ 6;
                sy = (k1y + 2.*k2y + 2.*k3y + k4y) ./ 6;

                x1 = x0 + sx .* dz;
                y1 = y0 + sy .* dz;
            end
        end

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
