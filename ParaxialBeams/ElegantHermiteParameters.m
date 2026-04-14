classdef ElegantHermiteParameters < GaussianParameters
    % ElegantHermiteParameters - Parameters for Elegant Hermite-Gaussian beams
    %
    % Mode indices:
    %   n  - horizontal mode index
    %   m  - vertical mode index
    %
    % Key quantity: complex beam parameter
    %   alpha(z) = i*k / (2*q(z))  where  q(z) = z + i*z_R
    %
    % Dynamic API (Phase 2, preferred for propagators):
    %   a   = params.alphaAtZ(z)   -- complex alpha at arbitrary z
    %   phi = params.phiPhase(z)   -- modal Gouy shift (n+m)*psi(z)
    %
    % Snapshot API (original, preserved):
    %   params.alpha     -- at stored zCoordinate
    %   params.PhiPhase  -- at stored zCoordinate

    properties
        n % Horizontal mode index
        m % Vertical mode index
    end

    properties (Dependent)
        alpha    % Complex beam parameter i*k/(2*q(z)) at stored zCoordinate
        PhiPhase % Modal Gouy phase (n+m)*psi(z) at stored zCoordinate
    end

    methods
        function obj = ElegantHermiteParameters(z, w0, lambda, n, m)
            % Constructor
            % z:      propagation distance (snapshot)
            % w0:     initial beam waist (m)
            % lambda: wavelength (m)
            % n, m:   mode indices (default 0)

            if nargin < 5
                m = 0;
            end
            if nargin < 4
                n = 0;
            end

            obj@GaussianParameters(z, w0, lambda);
            obj.n = n;
            obj.m = m;
        end

        % -----------------------------------------------------------------
        % Dynamic evaluation methods (Phase 2 API)
        % -----------------------------------------------------------------

        function a = alphaAtZ(obj, z)
            % alphaAtZ  Complex beam parameter alpha at arbitrary z.
            %
            % alpha(z) = i*k / (2*q(z))  where  q(z) = z + i*z_R
            %
            % This is the key parameter of the "elegant" variant: because
            % alpha is complex, the Hermite polynomial H_n(sqrt(alpha)*x)
            % evaluated with a complex argument produces amplitude AND phase
            % modulation simultaneously.
            a = obj.computeAlphaAtZ(z);
        end

        function phi = phiPhase(obj, z)
            % phiPhase  Modal Gouy phase shift (n+m)*psi(z) at arbitrary z.
            phi = (obj.n + obj.m) .* obj.gouyPhase(z);
        end

        % -----------------------------------------------------------------
        % Snapshot Dependent property getters (original API, preserved)
        % -----------------------------------------------------------------

        function a = get.alpha(obj)
            % Complex beam parameter at the stored zCoordinate.
            a = obj.computeAlpha();
        end

        function phi = get.PhiPhase(obj)
            % Modal Gouy phase shift at the stored zCoordinate.
            phi = (obj.n + obj.m) .* obj.GouyPhase;
        end
    end
end
