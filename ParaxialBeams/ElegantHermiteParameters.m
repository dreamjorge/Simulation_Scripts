classdef ElegantHermiteParameters < GaussianParameters
    % ElegantHermiteParameters - Parameters for Elegant Hermite-Gaussian beams
    
    properties
        n % Horizontal mode index
        m % Vertical mode index
    end
    
    properties (Dependent)
        alpha    % Complex beam parameter: i*k/(2*q(z)), q=z+i*zR
        PhiPhase % Mode Gouy phase: (n+m)*psi(z)
    end

    methods
        function obj = ElegantHermiteParameters(z, w0, lambda, n, m)
            % Constructor
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

        function a = get.alpha(obj)
            a = obj.computeAlpha();
        end

        function phi = get.PhiPhase(obj)
            phi = (obj.n + obj.m) .* obj.GouyPhase;
        end
    end
end
