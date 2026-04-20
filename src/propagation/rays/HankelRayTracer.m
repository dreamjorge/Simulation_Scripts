classdef HankelRayTracer < handle
    % HankelRayTracer — Eikonal ray tracing for Hankel (cylindrical) beams.
    %
    % ============================================================================
    % HANKEL BEAMS AND THEIR PHYSICS
    % ============================================================================
    %
    % Hankel beams are solutions of the paraxial wave equation in cylindrical
    % coordinates (r, theta, z).  Unlike Hermite-Gaussian beams (Cartesian),
    % they carry orbital angular momentum (OAM) and have a phase singularity
    % on axis.
    %
    % The complex amplitude of a Hankel-Laguerre beam is constructed as:
    %
    %   H_{lp}^{
  (1, 2)}(r, theta, z) = LG_{lp} +/- i * XLG_{lp}
    %
    % where LG_{lp} is the standard Laguerre-Gauss mode and XLG_{lp} is its
    % second independent solution (logarithmic + digamma series), analogous
    % to the Neumann function Y_l(kr) in Bessel theory.
    %
    %   H^{(1)} = LG + i*XLG   (outward propagating, diverging from axis)
    %   H^{(2)} = LG - i*XLG   (inward propagating, converging to axis)
    %
    % This superposition is the paraxial analog of the cylindrical Hankel
    % functions in free-space wave theory.  The XLG component carries the
    % radial energy flux information that determines whether the wave
    % radiates outward or converges inward.
    %
    % References:
    %   [Allen92]    L. Allen et al., "Orbital angular momentum of light
    %                and the transformation of Laguerre-Gaussian laser
    %                modes", Phys. Rev. A 45, 8185 (1992).
    %   [Kotlyar07]  V. V. Kotlyar et al., "Hankel-Bessel laser beams",
    %                J. Opt. Soc. Am. A 24, 1955-1966 (2007).
    %   [Ugalde20]   J. Ugalde, "Laguerre-Gauss beams and their Hankel
    %                decomposition", github.com/dreamjorge/LaguerreGaussBeams.
    %   [Papi23]     J. Papi, "Hankel beams", in Structured Light,
    %                Elsevier (2023).
    %
    % ============================================================================
    % EIKONAL RAY SLOPES: WHY (dphi/dz + k) MATTERS
    % ============================================================================
    %
    % Standard paraxial ray tracing computes slopes as:
    %
    %   sx = (1/k) * dphi/dx          (paraxial approximation)
    %
    % This assumes the longitudinal phase gradient is dominated by the
    % carrier: dphi/dz ~ -k, so the eikonal denominator reduces to k.
    %
    % For Hankel beams this approximation is WRONG.  Both H^(1) and H^(2)
    % share the same exp(-ikz) carrier and the same Gaussian envelope
    % exp(-r^2/w^2) * exp(ikr^2/2R).  Therefore the transverse gradient
    % dphi/dx is identical for both branches — the +-i*XLG modulation
    % affects only the AMPLITUDE, not the transverse phase, because XLG
    % enters multiplicatively inside the same carrier.
    %
    % The difference between H^(1) and H^(2) emerges in the LONGITUDINAL
    % gradient dphi/dz, where the z-evolution of the XLG component (via
    % the z-dependent Gouy phase shift, beam width w(z), and curvature
    % R(z)) creates a correction to the carrier-dominated dphi/dz.
    %
    % The full eikonal formula preserves this correction:
    %
    %   sx = (dphi/dx) / (dphi/dz + k)     (full eikonal)
    %   sy = (dphi/dy) / (dphi/dz + k)
    %
    % The +k removes the carrier contribution (since dphi/dz ~ -k +
    % envelope terms), isolating the envelope gradient in the denominator.
    % This denominator differs between H^(1) and H^(2), producing:
    %   - H^(1): positive radial slope (outward-diverging rays)
    %   - H^(2): negative radial slope (inward-converging rays)
    %
    % The paraxial formula divides by constant k ~ 10^7 rad/m, which
    % completely erases this small but physically essential difference.
    %
    % This is analogous to the approach in [Ugalde20], where the slope
    % ratio is computed as:
    %
    %   dz/dr = gz / gr,   with gz = gradient(phase_z)/dz + k
    %
    % and the ray step is dr = dz * gr / gz (full eikonal in cylindrical).
    %
    % References:
    %   [Born99]     M. Born & E. Wolf, Principles of Optics, 7th ed.,
    %                Cambridge (1999), Sec. 3.1.2 (eikonal equation).
    %   [Kravtsov90] Yu. A. Kravtsov & Yu. I. Orlov, Geometrical Optics
    %                of Inhomogeneous Media, Springer (1990), Ch. 2.
    %   [Ugalde20]   gradientrthz.m in LaguerreGaussBeams repository.
    %
    % ============================================================================
    % BRANCH SWITCHING: H^(1) <-> H^(2) AT THE AXIS
    % ============================================================================
    %
    % When an inward ray (H^(2)) reaches the optical axis, it must become
    % outward (H^(1)) and emerge on the OPPOSITE side.  This is the
    % cylindrical analog of a plane wave passing through a focus.
    %
    % The process has three parts:
    %
    %   1. CONVERGENCE — The full eikonal slopes drive H^(2) rays inward
    %      toward r = 0.  Without the correct eikonal, rays follow the
    %      Gaussian wavefront curvature outward and never reach the axis.
    %
    %   2. DETECTION — When the ray segment's minimum distance to the
    %      origin falls below a threshold, the crossing is detected.
    %      The threshold adapts to local geometry with a wavelength floor:
    %
    %        threshold = max( max(|x0|, |y0|) * 1e-3,  lambda )
    %
    %      The floor prevents the threshold from shrinking below detection
    %      range as the ray approaches the axis (race condition).
    %
    %   3. TRANSITION — The integration step is split at the crossing
    %      fraction t*.  After crossing:
    %        a) The Hankel type flips: H^(2) -> H^(1)
    %        b) The position reflects through origin: (x,y) -> (-x,-y)
    %        c) Integration continues with H^(1) slopes from the
    %           opposite side
    %
    % The position reflection matches the legacy cylindrical approach
    % in [Ugalde20] where ri < 0 triggers ri = abs(ri); xi = -xi.
    %
    % References:
    %   [Nienhuis96] G. Nienhuis, "Doppler effect induced by rotating
    %                lenses", Opt. Commun. 119, 271-285 (1996).
    %   [Courtial98] J. Courtial et al., "Gaussian beams with very high
    %                orbital angular momentum", Opt. Commun. 144 (1998).
    %   [Ugalde20]   getPropagateCylindricalRays in HankelLaguerre.m.
    %
    % ============================================================================
    % AXIS-CROSSING DETECTION: MINIMUM SEGMENT DISTANCE
    % ============================================================================
    %
    % For segment P(s) = P0 + s*(P1-P0) with s in [0,1]:
    %   s* = -P0 . (P1-P0) / |P1-P0|^2     (projection parameter)
    %   min_dist = |P0 + clamp(s*) * (P1-P0)|
    %
    % This replaces an earlier cross-product sign test, which produced
    % spurious flips for grazing trajectories that curved near the axis
    % without physically crossing it.
    %
    % References:
    %   [Eberly01]   D. Eberly, "Distance between point and line segment",
    %                Geometric Tools (2001).
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

            if nargin < 5 || isempty(method), method = 'RK4';
end if nargin < 4 || isempty(dzInternal), dzInternal = [];
end

    if isempty (zPlanes) bundleOut = bundleIn;
return;
end

    zPlanes = zPlanes( :).'; if any (diff(zPlanes) < 0)
                  error('zPlanes must be monotonically increasing.');
end

    % Current state from input bundle last slice x0 = bundleIn.x( :, :, end);
y0 = bundleIn.y( :, :, end);
z0 = bundleIn.z( :, :, end);
ht0 = bundleIn.ht( :, :, end);

zStart = z0(1, 1);
if abs (zPlanes(1) - zStart)
  > 1e-15 error('zPlanes(1) must match initial bundle z (%.6e).', zStart);
end

    % Output bundle starts at initial state
      only(fixed - plane samples) bundleOut = RayBundle(x0, y0, zStart);
bundleOut.ht = ht0;

            for
              kk = 2 : numel(zPlanes) zTarget = zPlanes(kk);
            if zTarget
              < zStart error('zPlanes must be increasing and >= initial z.');
                end

                % Internal step for this interval
                if isempty(dzInternal)
                    dzStep = max((zTarget - zStart) / 20, eps);
                else dzStep = dzInternal;
                end

                    % Propagate from current state to current
                          target using internal steps bTmp =
                    RayBundle(x0, y0, zStart);
                bTmp.ht = ht0;
                bTmp = HankelRayTracer.propagate(bTmp, beam, zTarget, dzStep,
                                                 method);

                x1 = bTmp.x( :, :, end);
                y1 = bTmp.y( :, :, end);
                z1 = bTmp.z( :, :, end);
                sx1 = bTmp.sx( :, :, end);
                sy1 = bTmp.sy( :, :, end);
                ht1 = bTmp.ht( :, :, end);

                bundleOut.addStep(x1, y1, z1, sx1, sy1, ht1);

                % advance state x0 = x1;
                y0 = y1;
                z0 = z1;
                ht0 = ht1;
                zStart = z0(1, 1);
            end
        end

        function bundle = propagate(bundle, beam, z_final, dz, method)
            % PROPAGATE — Integrate Hankel ray bundle with eikonal slopes
            % and axis-crossing detection.
            %
            % For HankelLaguerre beams, ray slopes are computed using the
            % full eikonal formula (see calculateSlopesEikonal) so that
            % H^(2) rays converge toward the axis and H^(1) rays diverge.
            %
            % When a ray segment passes within threshold distance of the
            % optical axis, the Hankel branch flips (H^(2) -> H^(1)) and
            % the position reflects through origin.  See class header for
            % the full detection and transition protocol.

            if nargin < 5, method = 'RK4';
            end

                z_current = bundle.z(1, 1, end);
            while
              z_current<z_final if z_current + dz> z_final dz =
                  z_final - z_current;
            end

                x0 = bundle.x( :, :, end);
            y0 = bundle.y( :, :, end);
            z0 = bundle.z( :, :, end);
            ht0 = bundle.ht( :, :, end);

            if strcmpi (method, 'Euler')
              [ sx, sy ] =
                  HankelRayTracer.calculateSlopes(beam, x0, y0, z0, ht0);
            x1 = x0 + sx.*dz;
            y1 = y0 + sy.*dz;
            else[k1x, k1y] =
                HankelRayTracer.calculateSlopes(beam, x0, y0, z0, ht0);
            [ k2x, k2y ] = HankelRayTracer.calculateSlopes(
                beam, x0 + k1x.*dz / 2, y0 + k1y.*dz / 2, z0 + dz / 2, ht0);
            [ k3x, k3y ] = HankelRayTracer.calculateSlopes(
                beam, x0 + k2x.*dz / 2, y0 + k2y.*dz / 2, z0 + dz / 2, ht0);
            [ k4x, k4y ] = HankelRayTracer.calculateSlopes(
                beam, x0 + k3x.*dz, y0 + k3y.*dz, z0 + dz, ht0);

            sx = (k1x + 2. * k2x + 2. * k3x + k4x)./ 6;
            sy = (k1y + 2. * k2y + 2. * k3y + k4y)./ 6;

            x1 = x0 + sx.*dz;
            y1 = y0 + sy.*dz;
            end

                z1 = z0 + dz;

            % Axis - crossing check — see class header[Eberly01] dx = x1 - x0;
            dy = y1 - y0;
            denominator = dx.^ 2 + dy.^ 2 + eps;

            % s * = -P₀·(P₁ - P₀) / | P₁ - P₀ |² t_param =
                        -(x0.*dx + y0.*dy)./ denominator;

            % Clamp projection to segment bounds t_clamped =
                max(0, min(1, t_param));

            % Nearest point on segment to origin nearest_x = x0 + t_clamped.*dx;
            nearest_y = y0 + t_clamped.*dy;

            % Minimum Euclidean distance min_dist =
                sqrt(nearest_x.^ 2 + nearest_y.^ 2);

            % Adaptive threshold with wavelength floor to prevent %
                the threshold from shrinking below detection range threshold =
                max(max(abs(x0), abs(y0)) * 1e-3, beam.Lambda);

            % Flip only if segment passes close enough to axis crossed =
                min_dist < threshold;

            ht1 = ht0;
            ht1(crossed &(ht0 == 2)) = 1;
            % flip H ^ (2) → H ^
                (1)

                        % Sub -
                    step crossing correction — see class header,
                % "BRANCH SWITCHING" section if any (crossed( :)) idxCross =
                    find(crossed);
                    for
                      ii = 1 : numel(idxCross) idx = idxCross(ii);

                    x0i = x0(idx);
                    y0i = y0(idx);
                    z0i = z0(idx);
                    ht0i = ht0(idx);

                    % Crossing fraction along this step(
                          already clamped[0, 1]) tCross = t_clamped(idx);
                    dz1 = dz * tCross;
                    dz2 = dz - dz1;

                    % Step 1 : up to crossing with original branch xMid = x0i;
                    yMid = y0i;
                    zMid = z0i;
                    if dz1
                      > 0 [xMid, yMid] = HankelRayTracer.singleStep(
                          beam, x0i, y0i, z0i, ht0i, dz1, method);
                    zMid = z0i + dz1;
                    end

                        % Step 2 : after crossing,
                        force H ^ (2)->H ^ (1)htMid = ht0i;
                    if ht0i
                      == 2 htMid = 1;
                    %
                        Reflect position through origin
                        : inward ray %
                          passed through axis and must emerge on the %
                          opposite side.Matches legacy behavior in %
                          HankelLaguerre.getPropagateCylindricalRays %
                          (xi = -xi; yi = -yi after H2->H1 flip).xMid = -xMid;
                    yMid = -yMid;
                    end

                        x1i = xMid;
                    y1i = yMid;
                    if dz2
                      > 0 [x1i, y1i] = HankelRayTracer.singleStep(
                          beam, xMid, yMid, zMid, htMid, dz2, method);
                    end

                        x1(idx) = x1i(1);
                    y1(idx) = y1i(1);
                    ht1(idx) = htMid;

                        % Effective slope stored for this global step
                        sx(idx) = (x1i(1) - x0i) / dz;
                        sy(idx) = (y1i(1) - y0i) / dz;
                        end end

                            bundle.addStep(x1, y1, z1, sx, sy, ht1);
                        z_current = z1(1, 1);
            end
        end


        function [sx, sy] = calculateSlopes(beam, x, y, z, ht)
            % CALCULATESLOPES — Ray slopes for Hankel beams, dispatching
            % to the appropriate gradient method by beam type.
            %
            % HankelLaguerre:
            %   Uses the full eikonal formula (calculateSlopesEikonal):
            %     sx = (dphi/dx) / (dphi/dz + k)
            %   This is required because the paraxial formula (1/k)*dphi/dx
            %   erases the XLG correction that differentiates H^(1) from
            %   H^(2).  See class documentation, "EIKONAL RAY SLOPES".
            %
            % Other beam types (HankelHermite, etc.):
            %   Falls back to the paraxial complex-field gradient via
            %   RayTracer.calculatePhaseGradientComplex.

            if isa(beam, 'HankelLaguerre')
                [sx, sy] = HankelRayTracer.calculateSlopesEikonal(beam, x, y, z, ht);
            return;
            end

                uniqueTypes = unique(ht( :));
            sx = zeros(size(x));
            sy = zeros(size(x));

            for
              t = 1 : numel(uniqueTypes) htype = uniqueTypes(t);
            mask = (ht == htype);

            tempBeam = HankelRayTracer.beamWithType(beam, htype);
            [ sx_part, sy_part ] =
                RayTracer.calculatePhaseGradientComplex(tempBeam, x, y, z);

            sx(mask) = sx_part(mask);
            sy(mask) = sy_part(mask);
            end end

                    end %
                methods(Static)

                    methods(Static, Access = private)

                        function[x1, y1] =
                singleStep(beam, x0, y0, z0, ht0, dz, method) %
                    SINGLESTEP — Integrate one ray over a scalar dz.%
                    Used by crossing sub
                - step correction to split H2->H1 transition.

                  if strcmpi (method, 'Euler')[sx, sy] =
                    HankelRayTracer.calculateSlopes(beam, x0, y0, z0, ht0);
            x1 = x0 + sx.*dz;
            y1 = y0 + sy.*dz;
            else[k1x, k1y] =
                HankelRayTracer.calculateSlopes(beam, x0, y0, z0, ht0);
            [ k2x, k2y ] = HankelRayTracer.calculateSlopes(
                beam, x0 + k1x.*dz / 2, y0 + k1y.*dz / 2, z0 + dz / 2, ht0);
            [ k3x, k3y ] = HankelRayTracer.calculateSlopes(
                beam, x0 + k2x.*dz / 2, y0 + k2y.*dz / 2, z0 + dz / 2, ht0);
            [ k4x, k4y ] = HankelRayTracer.calculateSlopes(
                beam, x0 + k3x.*dz, y0 + k3y.*dz, z0 + dz, ht0);

            sx = (k1x + 2. * k2x + 2. * k3x + k4x)./ 6;
            sy = (k1y + 2. * k2y + 2. * k3y + k4y)./ 6;

            x1 = x0 + sx.*dz;
            y1 = y0 + sy.*dz;
            end
        end

        function [sx, sy] = calculateSlopesEikonal(beam, x, y, z, ht)
            % CALCULATESLOPESEIKONAL — Full eikonal slopes for Hankel beams.
            %
            % Computes ray slopes using the exact eikonal relation instead
            % of the paraxial approximation:
            %
            %   Paraxial:  sx = (1/k) * dphi/dx              (WRONG for Hankel)
            %   Eikonal:   sx = (dphi/dx) / (dphi/dz + k)    (CORRECT)
            %
            % WHY THE PARAXIAL FORMULA FAILS:
            %
            %   The Hankel field is H = LG +/- i*XLG, where both LG and
            %   XLG share the same Gaussian carrier exp(-ikz)*exp(ikr^2/2R).
            %   The transverse gradient dphi/dx is dominated by this shared
            %   carrier — the +/-i*XLG modulates the amplitude, not the
            %   transverse phase.  Dividing by the constant k ~ 10^7 rad/m
            %   produces slopes that are identical for H^(1) and H^(2) to
            %   machine precision (~1e-15 relative difference).
            %
            %   The LONGITUDINAL gradient dphi/dz differs between H^(1)
            %   and H^(2) because the XLG component has z-dependent terms
            %   (beam width w(z), curvature R(z), Gouy phase) that evolve
            %   differently from LG.  The ratio dphi_perp / (dphi/dz + k)
            %   captures this difference, producing:
            %     H^(1): positive radial component (diverging)
            %     H^(2): negative radial component (converging)
            %
            % NUMERICAL METHOD:
            %
            %   All three phase gradients (dphi/dx, dphi/dy, dphi/dz) are
            %   computed via the complex-field identity:
            %
            %     dphi/dq = Im{ u* * du/dq } / |u|^2      [Takeda84]
            %
            %   with du/dq approximated by central difference.  The "+k"
            %   in the denominator removes the carrier exp(-ikz) whose
            %   gradient is -k (our sign convention).
            %
            %   For vortex beams (l != 0), the transverse derivatives are
            %   computed in polar coordinates (r, theta) and transformed
            %   to Cartesian, avoiding branch-cut artifacts at theta = 0.
            %
            % RELATION TO LEGACY CYLINDRICAL APPROACH:
            %
            %   The legacy code (gradientrthz.m in [Ugalde20]) computes:
            %
            %     gz = gradient(phase_z)/dz + k
            %     gr = gradient(phase_r)/dr
            %     dz/dr = gz / gr   =>   dr/dz = gr / gz
            %
            %   This is the same eikonal in cylindrical coordinates.  Our
            %   implementation generalizes it to 2D Cartesian (x,y) with
            %   RK4 integration, supporting both vortex and non-vortex
            %   modes uniformly.
            %
            % References:
            %   [Born99]     M. Born & E. Wolf, Principles of Optics,
            %                7th ed., Cambridge (1999), Sec. 3.1.2.
            %   [Takeda84]   M. Takeda et al., "Fourier-transform method
            %                of fringe-pattern analysis", J. Opt. Soc.
            %                Am. 72, 156-160 (1982).
            %   [Kravtsov90] Yu. A. Kravtsov & Yu. I. Orlov, Geometrical
            %                Optics of Inhomogeneous Media, Springer (1990).
            %   [Ugalde20]   gradientrthz.m in LaguerreGaussBeams.

            epsilon = 1e-12;
            w0 = beam.InitialWaist;
            lambda = beam.Lambda;
            k = beam.k;

            uniqueTypes = unique(ht( :));
            sx = zeros(size(x));
            sy = zeros(size(x));

            for
              t = 1 : numel(uniqueTypes) htype = uniqueTypes(t);
            mask = (ht == htype);

            tempBeam = HankelRayTracer.beamWithType(beam, htype);

            delta_mat = RayTracer.resolveDelta(x, y, w0, lambda);
            delta = max(delta_mat( :));
            % Compute the ENVELOPE phase gradient in z by stripping %
                    the carrier exp(-ikz) before differentiating.%
                    % Direct central difference of dphidz
                + k suffers from % catastrophic cancellation : dphidz ~-k +
                O(1 / zr),
                so % dphidz + k ~O(1 / zr) is lost in the ~k noise.Any % finite
                    - difference error in estimating -
                    k propagates % directly into gz.%
                        % Instead : multiply u by exp(+ikz) to get the envelope,
                % then differentiate the envelope phase.The result IS
                    % gz = dphidz + k without cancellation.zr_beam = pi * w0 ^
                                                                     2 / lambda;
            dz_z = max(lambda, max(abs(z( :))) * 1e-4);
            dz_z = max(dz_z, zr_beam * 1e-6);

            u0 = tempBeam.opticalField(x, y, z);
            u_zp = tempBeam.opticalField(x, y, z + dz_z);
            u_zm = tempBeam.opticalField(x, y, z - dz_z);

            % Strip carrier
                : u_env(z) = u(z) * exp(+ikz) u0_env = u0.*exp(1i * k * z);
            u_zp_env = u_zp.*exp(1i * k * (z + dz_z));
            u_zm_env = u_zm.*exp(1i * k * (z - dz_z));

            u0c_env = conj(u0_env);
            abs_u0_sq = real(u0c_env.*u0_env) + epsilon;

            % gz = d(phase_envelope) / dz =
                       dphidz + k(carrier removed) gz =
                           imag(u0c_env.*(u_zp_env - u_zm_env)./ (2 * dz_z))./
                           abs_u0_sq;
            gz(abs(gz) < epsilon) = epsilon;

            u0c = conj(u0);
delta);
[ x_tp, y_tp ] = pol2cart(TH + delta_theta, R);
[ x_tm, y_tm ] = pol2cart(TH - delta_theta, R);

u_rp = tempBeam.opticalField(x_rp, y_rp, z);
u_rm = tempBeam.opticalField(x_rm, y_rm, z);
u_tp = tempBeam.opticalField(x_tp, y_tp, z);
u_tm = tempBeam.opticalField(x_tm, y_tm, z);

dphidr = imag(u0c.*(u_rp - u_rm)./ (2 * delta))./ abs_u0_sq;
dphidt = imag(u0c.*(u_tp - u_tm)./ (2.*delta_theta))./ abs_u0_sq;

R_reg = R + eps;
R_sq = R.^ 2 + eps;
sx_part = (dphidr.*x./ R_reg - dphidt.*y./ R_sq)./ gz;
sy_part = (dphidr.*y./ R_reg + dphidt.*x./ R_sq)./ gz;
else u_xp = tempBeam.opticalField(x + delta, y, z);
u_xm = tempBeam.opticalField(x - delta, y, z);
u_yp = tempBeam.opticalField(x, y + delta, z);
u_ym = tempBeam.opticalField(x, y - delta, z);

dphidx = imag(u0c.*(u_xp - u_xm)./ (2 * delta))./ abs_u0_sq;
dphidy = imag(u0c.*(u_yp - u_ym)./ (2 * delta))./ abs_u0_sq;

sx_part = dphidx./ gz;
sy_part = dphidy./ gz;
end

    sx(mask) = sx_part(mask);
sy(mask) = sy_part(mask);
end end

    function newBeam = beamWithType(beam, htype) %
                       BEAMWITHTYPE — Create
                           a beam copy with a specific Hankel branch.%
                       % Hankel beams carry a type index : 1 = H ^ (1)(outward),
             % 2 = H ^
                       (2)(inward).This method fabricates the %
                               correct beam variant without re -
                           computing the field formula.% %
                               For HankelLaguerre
    : the Hankel type determines whether %
      the radial Bessel function uses H ^
                       (1) or
                   H ^ (2),
             which % have different asymptotic behavior at infinity.% %
                 For Hermite beams
    : no branch
      switching(Cartesian,
                % no singular axis) — returns original beam unchanged.

      if isa (beam, 'HankelLaguerre') newBeam = HankelLaguerre(
                 beam.InitialWaist, beam.Lambda, ... beam.l, beam.p, htype);
elseif isa(beam, 'HankelHermite') newBeam = HankelHermite(beam.InitialWaist,
                                                          beam.Lambda,
                                                          ... beam.n, beam.m,
                                                          htype);
else newBeam = beam;
end end

        end %
    methods(Static, Access = private)

        end
    % classdef
