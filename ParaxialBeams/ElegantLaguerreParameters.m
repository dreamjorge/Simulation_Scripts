classdef ElegantLaguerreParameters < GaussianParameters
    % ElegantLaguerreParameters - Parameters for Elegant Laguerre-Gaussian beams
    %
    % These beams use a complex argument in the Laguerre polynomial, scaled
    % by the complex beam parameter alpha instead of the real waist w.
    %
    % Mode indices:
    %   l  - topological charge (azimuthal index, integer, may be negative)
    %   p  - radial order (non-negative integer)
    %
    % Key quantity: complex beam parameter
    %   alpha(z) = i*k / (2*q(z))  where  q(z) = z + i*z_R
    %
    % Dynamic API (Phase 2, preferred for propagators):
    %   a   = params.alphaAtZ(z)   -- complex alpha at arbitrary z
    %   phi = params.phiPhase(z)   -- modal Gouy shift (|l|+2p)*psi(z)
    %
    % Snapshot API (original, preserved):
    %   params.alpha     -- at stored zCoordinate
    %   params.PhiPhase  -- at stored zCoordinate

    properties
        l % Topological charge (azimuthal index)
        p % Radial index
    end

    properties (Dependent)
        alpha    % Complex beam parameter i*k/(2*q(z)) at stored zCoordinate
        PhiPhase % Modal Gouy phase (|l|+2p)*psi(z) at stored zCoordinate
    end

    methods
        function obj = ElegantLaguerreParameters(z, w0, lambda, l, p)
            % Constructor
            % z:      propagation distance (snapshot)
            % w0:     initial beam waist (m)
            % lambda: wavelength (m)
            % l:      topological charge (default 0)
            % p:      radial order (default 0)

            if nargin < 5
                p = 0;
            end
            if nargin < 4
                l = 0;
            end

            obj@GaussianParameters(z, w0, lambda);
            obj.l = l;
            obj.p = p;
        end

        % -----------------------------------------------------------------
        % Dynamic evaluation methods (Phase 2 API)
        % -----------------------------------------------------------------

        function a = alphaAtZ(obj, z)
            % alphaAtZ  Complex beam parameter alpha at arbitrary z.
            %
            % alpha(z) = i*k / (2*q(z))  where  q(z) = z + i*z_R
            %
            % This is the key parameter of the "elegant" variant: using alpha
            % instead of the real waist w changes both the radial amplitude
            % (sqrt(alpha)*r)^|l| and the polynomial argument alpha*r^2,
            % producing the characteristic elegant amplitude-phase coupling.
            a = obj.computeAlphaAtZ(z);
        end

        function phi = phiPhase(obj, z)
            % phiPhase  Modal Gouy phase shift (|l|+2p)*psi(z) at arbitrary z.
            phi = (abs(obj.l) + 2*obj.p) .* obj.gouyPhase(z);
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
            phi = (abs(obj.l) + 2*obj.p) .* obj.GouyPhase;
        end
    end
end
