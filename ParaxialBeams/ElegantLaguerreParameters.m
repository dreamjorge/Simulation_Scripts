classdef ElegantLaguerreParameters < GaussianParameters
    % ElegantLaguerreParameters - Parameters for Elegant Laguerre-Gaussian beams
    % These beams use a complex argument in the Laguerre polynomial.
    
    properties
        l % Topological charge (angular index)
        p % Radial index
    end
    
    properties (Dependent)
        alpha % Complex parameter for the Laguerre polynomial argument
    end
    
    methods
        function obj = ElegantLaguerreParameters(z, w0, lambda, l, p)
            % Constructor
            % z: axial distance
            % w0: initial waist
            % lambda: wavelength
            % l: topological charge
            % p: radial index
            
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
            % alpha = i*k / (2*q(z))
            % q(z) = z + i*zr
            k_val = obj.k;
            zr = obj.RayleighDistance;
            z_val = obj.zCoordinate;
            q = z_val + 1i * zr;
            a = 1i * k_val / (2 * q);
        end
    end
end
