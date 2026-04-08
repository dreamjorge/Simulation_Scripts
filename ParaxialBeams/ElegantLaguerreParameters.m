classdef ElegantLaguerreParameters < GaussianParameters
    % ElegantLaguerreParameters - Parameters for Elegant Laguerre-Gaussian beams
    % These beams use a complex argument in the Laguerre polynomial.
    
    properties
        l % Topological charge (angular index)
        p % Radial index
    end
    
    properties (Dependent)
        alpha    % Complex beam parameter: i*k/(2*q(z)), q=z+i*zR
        PhiPhase % Mode Gouy phase: (|l|+2p)*psi(z)
    end

    methods
        function obj = ElegantLaguerreParameters(z, w0, lambda, l, p)
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

        function a = get.alpha(obj)
            a = obj.computeAlpha();
        end

        function phi = get.PhiPhase(obj)
            phi = (abs(obj.l) + 2*obj.p) .* obj.GouyPhase;
        end
    end
end
