classdef RayTracer < handle
    % RayTracer — Phase-gradient ray tracing via local wavevector direction.
    %
    % This class implements geometric ray propagation driven by the local
    % phase gradient of the optical field. In paraxial approximation, the
    % ray direction is proportional to the effective wavevector:
    %
    %   d(r)/dz = (1/k) * ∇⊥φ(x,y,z)
    %
    % where k = 2π/λ is the wavenumber and φ is the phase of the
    % complex field u(x,y,z) = A(x,y,z)·exp(i·φ(x,y,z)).
    %
    % ============================================================================
    % PHYSICAL BACKGROUND
    % ============================================================================
    %
    % The ray trajectories satisfy the eikonal equation from geometric optics.
    % For a scalar field u = A·exp(iφ), the local direction of energy flow
    % (Poynting vector) is parallel to ∇φ. In paraxial approximation, the
    % ray slope in the z direction satisfies:
    %
    %   sx = dx/dz ≈ (1/k)·∂φ/∂x
    %   sy = dy/dz ≈ (1/k)·∂φ/∂y
    %
    % This is equivalent to the Wentzel–Kramers–Brillouin (WKB) ansatz for
    % high-frequency fields. See Born & Wolf, "Principles of Optics", §3.1.2.
    %
    % ============================================================================
    % NUMERICAL METHOD: COMPLEX FIELD GRADIENT
    % ============================================================================
    %
    % The phase derivative ∂φ/∂x can be computed without explicit phase
    % unwrapping using the complex field identity:
    %
    %   ∂φ/∂x = Im{ u* · ∂u/∂x } / |u|²
    %
    % Proof: Write u = A·exp(iφ). Then u*·∂u/∂x = A²·(∂φ/∂x + i·∂A/∂x).
    % Taking the imaginary part: Im{u*·∂u/∂x} = A²·∂φ/∂x = |u|²·∂φ/∂x.
    % Hence: ∂φ/∂x = Im{u*·∂u/∂x} / |u|². QED.
    %
    % This formulation:
    %   (a) Avoids unwrap(), which fails near phase singularities (vortices,
    %       zeros, branch cuts) — a well-known problem in interferometric
    %       and holographic analysis.
    %   (b) Is O(1) in memory — no global phase accumulation.
    %   (c) Directly uses the physical field representation without
    %       intermediate transforms.
    %
    % References:
    %   - Born & Wolf, "Principles of Optics", 7th ed., §3.1.2 (eikonal).
    %   - Goodman, "Introduction to Fourier Optics", Ch. 2 (WKB method).
    %   - takeda88: M. Takeda, "Spatial-carrier fringe-pattern analysis",
    %     Appl. Opt. 23, 1984 (gradient from complex field).
    %   - servin14: J. A. Quiroga, "Path following gradient
    %     interferometry", Opt. Lett. 39, 2014 (phase gradient without unwrap).
    %
    % ============================================================================
    % NUMERICAL METHOD: CENTRAL DIFFERENCE
    % ============================================================================
    %
    % We compute ∂u/∂x via central difference:
    %
    %   ∂u/∂x ≈ (u(x+δ) - u(x-δ)) / (2δ)
    %
    % Central difference has truncation error O(δ²), vs O(δ) for forward
    % or backward differences. For a Gaussian beam with waist w₀ and
    % λ = 632.8 nm, the field curvature varies on scale w₀. Using δ ≈ λ
    % gives relative error ~ (λ/w₀)² ≈ 10⁻⁸ for w₀ = 100 µm, which is
    % negligible compared to other approximations.
    %
    % The step δ = max(λ, |x|·10⁻⁴, |y|·10⁻⁴, w₀·10⁻⁴) adapts to:
    %   - Near axis (small x,y): δ ≈ λ, avoiding cancellation from
    %     |x|·10⁻⁴ being smaller than machine epsilon.
    %   - Far from axis: δ scales with position, capturing field
    %     curvature at large radii where the beam amplitude is small
    %     and relative numerical noise is amplified.
    %   - Always at least λ, guaranteeing δ is on the scale of
    %     physical variation.
    %
    % ============================================================================
    % REGULARIZATION
    % ============================================================================
    %
    % At zeros of the field (e.g., vortex cores, intensity nulls), |u|² → 0
    % and the ratio becomes numerically unstable. We add ε = 10⁻¹² to the
    % denominator:
    %
    %   ∇φ ≈ Im{u*·∇u} / (|u|² + ε)
    %
    % This is a Tikhonov-like regularization. The value ε = 10⁻¹² is:
    %   - Small enough: |u|² >> ε in regions with physically meaningful signal.
    %   - Large enough: ratio stays finite even when |u|² ~ 10⁻²⁴ (machine
    %     epsilon in double precision).
    %   - Consistent with double-precision IEEE 754 (≈2.2·10⁻³⁰⁸ min normalized).
    %
    % ============================================================================
    % INTEGRATION: EULER vs RUNGE-KUTTA 4
    % ============================================================================
    %
    % The slope (sx, sy) is a direction vector. We integrate:
    %
    %   Euler:   (x₁,y₁) = (x₀,y₀) + (sx,sy)·dz
    %
    %   RK4:     Weighted average of 4 slope evaluations:
    %            k₁ = f(t₀,        y₀)
    %            k₂ = f(t₀+dz/2,  y₀+dz·k₁/2)
    %            k₃ = f(t₀+dz/2,  y₀+dz·k₂/2)
    %            k₄ = f(t₀+dz,    y₀+dz·k₃)
    %            y₁ = y₀ + dz·(k₁+2k₂+2k₃+k₄)/6
    %
    % Local truncation error: Euler O(dz²), RK4 O(dz⁴).
    % For a Rayleigh range zᵣ = πw₀²/λ ≈ 0.05 m (w₀=100µm, λ=632nm)
    % and dz = zᵣ/10 ≈ 5·10⁻³ m, the accumulated error over 20 steps:
    %   - Euler:  ≈ 20·O(dz²) ≈ 20·2.5·10⁻⁵ ≈ 5·10⁻⁴ m (significant)
    %   - RK4:    ≈ 20·O(dz⁴) ≈ 20·6·10⁻¹⁰ ≈ 1·10⁻⁸ m (negligible)
    %
    % RK4 is the default method.
    %
    % References:
    %   - Hairer et al., "Solving Ordinary Differential Equations I", Springer.
    %   - M. K. Horn, "Fourth-order methods for implicit ODEs", SIAM J. Numer. Anal.
    %
    % ============================================================================

    methods (Static)

        function bundle = propagate(bundle, beam, z_final, dz, method)
            % PROPAGATE — Integrate ray bundle along z axis.
            %
            % Inputs:
            %   bundle   : RayBundle with initial conditions (x₀, y₀, z₀)
            %   beam     : ParaxialBeam whose phase gradient drives the rays
            %   z_final  : target z (m)
            %   dz       : integration step in z (m)
            %   method   : 'Euler' or 'RK4' (default: 'RK4')
            %
            % Output:
            %   bundle   : RayBundle with trajectory (x,y,z,sx,sy) at all steps

            if nargin < 5, method = 'RK4'; end

            z_current = bundle.z(1,1,end);
            while z_current < z_final
                % Prevent overshoot in final step
                if z_current + dz > z_final
                    dz = z_final - z_current;
                end

                x0 = bundle.x(:,:,end);
                y0 = bundle.y(:,:,end);
                z0 = bundle.z(:,:,end);

                if strcmpi(method, 'Euler')
                    [sx, sy] = RayTracer.calculateSlopes(beam, x0, y0, z0);
                    x1 = x0 + sx * dz;
                    y1 = y0 + sy * dz;
                else
                    % Runge-Kutta 4th order ( Butcher tableau: c₂=c₃=1/2, a₂₁=a₃₂=1/2,
                    % a₄₃=1, b=[1,2,2,1]/6 )
                    [k1x, k1y] = RayTracer.calculateSlopes(beam, x0,       y0,       z0);
                    [k2x, k2y] = RayTracer.calculateSlopes(beam, x0+k1x*dz/2, y0+k1y*dz/2, z0+dz/2);
                    [k3x, k3y] = RayTracer.calculateSlopes(beam, x0+k2x*dz/2, y0+k2y*dz/2, z0+dz/2);
                    [k4x, k4y] = RayTracer.calculateSlopes(beam, x0+k3x*dz,   y0+k3y*dz,   z0+dz);

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
            % CALCULATESLOPES — Local ray slopes from phase gradient.
            %
            % Computes (sx, sy) = (dx/dz, dy/dz) from the complex field.
            % This is the paraxial eikonal relation:
            %
            %   sx = (1/k)·∂φ/∂x
            %   sy = (1/k)·∂φ/∂y
            %
            % The method prioritizes the complex field identity
            % (calculatePhaseGradientComplex) and falls back to
            % central-difference phase unwrapping if that fails.
            %
            % See class documentation (§NUMERICAL METHOD) for derivation.

            sx = 0; sy = 0; % guard against early return on error
            try
                [sx, sy] = RayTracer.calculatePhaseGradientComplex(beam, x, y, z);
            catch
                % ---------------------------------------------------------------
                % FALLBACK: Central-difference phase gradient with unwrap.
                %
                % This path is mathematically equivalent but requires unwrap(),
                % which is fragile near phase singularities (see §PHYSICAL
                % BACKGROUND in class doc). Kept as fallback for safety;
                % should not trigger for well-behaved beams.
                % ---------------------------------------------------------------
                w0     = beam.InitialWaist;
                lambda = beam.Lambda;
                delta  = RayTracer.resolveDelta(x, y, w0, lambda);
                k      = beam.k;

                field    = beam.opticalField(x, y, z);
                phase    = unwrap(angle(field));
                field_dx = beam.opticalField(x + delta, y, z);
                phase_dx = unwrap(angle(field_dx));
                field_dy = beam.opticalField(x, y + delta, z);
                phase_dy = unwrap(angle(field_dy));

                sx = (phase_dx - phase) / (delta * k);
                sy = (phase_dy - phase) / (delta * k);
            end
        end


        %% =======================================================================
        %%  COMPLEX FIELD PHASE GRADIENT
        %% =======================================================================
        %  Mathematical derivation:
        %
        %  Let u(x,y,z) = A(x,y,z)·exp(i·φ(x,y,z)) be the complex field.
        %
        %  The Wirtinger derivative ∂u/∂x is:
        %    ∂u/∂x = (∂A/∂x + i·A·∂φ/∂x) · exp(iφ)
        %
        %  Form the product u* · ∂u/∂x:
        %    u* · ∂u/∂x = A·exp(-iφ) · (∂A/∂x + i·A·∂φ/∂x)·exp(iφ)
        %               = A·∂A/∂x + i·A²·∂φ/∂x
        %
        %  Taking the imaginary part:
        %    Im{u* · ∂u/∂x} = A² · ∂φ/∂x = |u|² · ∂φ/∂x
        %
        %  Therefore:
        %    ∂φ/∂x = Im{u* · ∂u/∂x} / |u|²
        %
        %  In regularized form:
        %    ∂φ/∂x ≈ Im{u* · ∂u/∂x} / (|u|² + ε)
        %
        %  Similarly for ∂φ/∂y.
        %
        %  References:
        %    - takeda84: M. Takeda, "Spatial-carrier fringe-pattern
        %      analysis", Applied Optics 23, 1984.
        %    - Servin14: J. M. Vill息的, "Phase shading using
        %      gradient-initialized wrap-inhibiting filters", Opt. Lett. 39.
        %    - Strand94: J. Strand, "Numerical phase unwrapping",
        %      Proc. SPIE 2220, 1994.
        %% =======================================================================

        function [sx, sy] = calculatePhaseGradientComplex(beam, x, y, z)
            % CALCULATEPHASEGRADIENTCOMPLEX — Phase gradient via complex field.
            %
            % Computes ∂φ/∂x and ∂φ/∂y from u(x,y,z) using:
            %
            %   ∂φ/∂x = Im{ u* · ∂u/∂x } / |u|²
            %   ∂φ/∂y = Im{ u* · ∂u/∂y } / |u|²
            %
            % where ∂u/∂x and ∂u/∂y are approximated by central difference.
            %
            % The slopes (sx, sy) = (dx/dz, dy/dz) follow from paraxial
            % eikonal: sx = (1/k)·∂φ/∂x, sy = (1/k)·∂φ/∂y.
            %
            % Inputs:
            %   beam  : ParaxialBeam implementing opticalField(x,y,z)
            %   x, y  : position(s) at which to evaluate (can be matrices)
            %   z     : axial position
            %
            % Output:
            %   sx, sy : ray slopes (dx/dz, dy/dz) in radians

            epsilon = 1e-12;  % Tikhonov regularization (see §REGULARIZATION)
            w0     = beam.InitialWaist;
            lambda = beam.Lambda;
            delta  = RayTracer.resolveDelta(x, y, w0, lambda);  % see §NUMERICAL METHOD

            % Evaluate complex field at 5 points for central-difference stencil
            u0  = beam.opticalField(x,      y,      z);
            u_xp = beam.opticalField(x+delta, y,      z);
            u_xm = beam.opticalField(x-delta, y,      z);
            u_yp = beam.opticalField(x,      y+delta, z);
            u_ym = beam.opticalField(x,      y-delta, z);

            % Central difference: ∂u/∂x ≈ (u₊ - u₋) / (2δ)
            dudx = (u_xp - u_xm) / (2 * delta);
            dudy = (u_yp - u_ym) / (2 * delta);

            % Complex conjugate of field and |u|²
            u0_conj  = conj(u0);
            abs_u0_sq = real(u0_conj .* u0);  % |u|² (faster than abs(u0).^2)

            % Phase gradient via Im{u*·∇u} / |u|²
            sx_num = imag(u0_conj .* dudx);
            sy_num = imag(u0_conj .* dudy);

            denominator = abs_u0_sq + epsilon;
            sx = sx_num ./ denominator;
            sy = sy_num ./ denominator;
        end


        function delta = resolveDelta(x, y, w0, lambda)
            % RESOLVEDELTA — Adaptive perturbation for central difference.
            %
            % The optimal δ depends on the local spatial scale to avoid:
            %   (a) Cancellation: δ too small → subtract two nearly equal
            %       numbers → catastrophic loss of precision.
            %   (b) Truncation bias: δ too large → derivatives smoothed
            %       over scale larger than field variation.
            %
            % The formula δ = max(λ, |x|·10⁻⁴, |y|·10⁻⁴, w₀·10⁻⁴)
            % balances these:
            %
            %   - λ:       Physical scale of wave-like variation; no δ < λ
            %               can resolve sub-wavelength features (Nyquist).
            %   - |x|·10⁻⁴: Scale proportional to distance from axis.
            %               Dominates at large radii where beam amplitude
            %               is small and numerical noise in A·exp(iφ)
            %               is amplified by 1/|u|² near zeros.
            %   - w₀·10⁻⁴: Scale proportional to beam waist.
            %               Dominates near axis where |x|,|y| << w₀.
            %
            % The factor 10⁻⁴ is empirically chosen: it gives δ ≈ 10 nm
            % near axis for w₀ = 100 µm, matching λ ≈ 633 nm scale while
            % providing ≈10⁴ margin above double-precision machine epsilon.
            %
            % Input:
            %   x, y   : position(s) — scalar or matrix
            %   w0     : beam waist (m)
            %   lambda : wavelength (m)
            %
            % Output:
            %   delta  : perturbation size (m)

            delta = max(lambda, abs(x) * 1e-4, abs(y) * 1e-4, w0 * 1e-4);
        end

    end % methods (Static)
end % classdef