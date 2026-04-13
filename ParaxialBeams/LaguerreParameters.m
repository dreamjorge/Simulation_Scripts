classdef LaguerreParameters < GaussianParameters
    % LaguerreParameters - Parameters for Laguerre-Gaussian beams
    % Inherits from GaussianParameters and adds Laguerre-specific properties.
    %
    % Mode indices:
    %   l  - topological charge (azimuthal index, integer, may be negative)
    %   p  - radial order (non-negative integer)
    %
    % Dynamic API (Phase 2, preferred for propagators):
    %   phi = params.phiPhase(z)       -- modal Gouy shift (|l|+2p)*psi(z)
    %   w   = params.laguerreWaist(z)  -- spot size w(z)*sqrt(2p+|l|+1)
    %
    % Snapshot API (original, preserved):
    %   params.PhiPhase      -- at stored zCoordinate
    %   params.LaguerreWaist -- at stored zCoordinate

    properties
        l % Topological charge (azimuthal index)
        p % Radial order
    end

    properties (Dependent)
        LaguerreWaist  % Spot size of Laguerre-Gaussian beam at stored z
        PhiPhase       % Modal Gouy phase shift at stored z
    end

    methods
        function obj = LaguerreParameters(z, w0, lambda, l, p)
            % Constructor
            % z:      propagation distance (snapshot)
            % w0:     initial beam waist (m)
            % lambda: wavelength (m)
            % l:      topological charge (default 0)
            % p:      radial order (default 0)

            if nargin < 3
                error('Usage: params = LaguerreParameters(z, w0, lambda, l, p)');
            end

            obj@GaussianParameters(z, w0, lambda);

            if nargin >= 5
                obj.l = l;
                obj.p = p;
            else
                obj.l = 0;
                obj.p = 0;
            end
        end

        % -----------------------------------------------------------------
        % Dynamic evaluation methods (Phase 2 API)
        % -----------------------------------------------------------------

        function phi = phiPhase(obj, z)
            % phiPhase  Modal Gouy phase shift (|l|+2p)*psi(z) at arbitrary z.
            %
            % This is the additional Gouy phase accumulated by mode LG_{lp}
            % relative to the fundamental Gaussian mode. At z = 0 it is zero.
            % The sign of l does not affect the Gouy shift (abs value is used).
            phi = (abs(obj.l) + 2*obj.p) .* obj.gouyPhase(z);
        end

        function wL = laguerreWaist(obj, z)
            % laguerreWaist  LG spot size w(z)*sqrt(2p+|l|+1) at arbitrary z.
            wL = obj.waist(z) .* sqrt(2*obj.p + abs(obj.l) + 1);
        end

        % -----------------------------------------------------------------
        % Snapshot Dependent property getters (original API, preserved)
        % -----------------------------------------------------------------

        function phi = get.PhiPhase(obj)
            % Modal Gouy phase shift (|l|+2p)*psi at stored zCoordinate.
            phi = (abs(obj.l) + 2*obj.p) .* obj.GouyPhase;
        end

        function wL = get.LaguerreWaist(obj)
            % Spot size of LG beam: w(z)*sqrt(2p+|l|+1).
            wL = obj.Waist .* sqrt(2*obj.p + abs(obj.l) + 1);
        end
    end

    methods (Static)
        function wL = getWaist(z, w0, zr, l, p)
            % Static: Laguerre-Gaussian waist at z.
            w  = w0 * sqrt(1 + (z/zr).^2);
            wL = w * sqrt(2*p + abs(l) + 1);
        end
    end
end
