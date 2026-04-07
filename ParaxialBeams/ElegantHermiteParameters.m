classdef ElegantHermiteParameters < GaussianParameters
    % ElegantHermiteParameters - Parameters for Elegant Hermite-Gaussian beams
    
    properties
        n % Horizontal mode index
        m % Vertical mode index
    end
    
    properties (Dependent)
        alpha % Complex parameter
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
            k_val = obj.k;
            zr = obj.RayleighDistance;
            z_val = obj.zCoordinate;
            q = z_val + 1i * zr;
            a = 1i * k_val / (2 * q);
        end
    end
end
