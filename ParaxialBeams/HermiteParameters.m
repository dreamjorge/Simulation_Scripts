classdef HermiteParameters < GaussianParameters
    % HermiteParameters - Parameters for Hermite-Gaussian beams
    % Inherits from GaussianParameters and adds Hermite-specific properties
    
    properties
        n   % Order in x
        m   % Order in y
    end
    
    properties (Dependent)
        HermiteWaistX
        HermiteWaistY
        HermiteWaist
        PhiPhase        % Total phase shift including Hermite orders
    end
    
    methods
        function obj = HermiteParameters(z, w0, lambda, n, m)
            % Constructor
            % z: propagation distance
            % w0: initial waist
            % lambda: wavelength
            % n, m: Hermite orders
            
            if nargin < 3
                error('Usage: params = HermiteParameters(z, w0, lambda, n, m)');
            end
            
            obj@GaussianParameters(z, w0, lambda);
            
            if nargin >= 5
                obj.n = n;
                obj.m = m;
            else
                obj.n = 0;
                obj.m = 0;
            end
        end
        
        %% Dependent property getters
        function phi = get.PhiPhase(obj)
            % Total phase shift: (n + m + 1) * psi(z)
            % Note: the +1 is usually handled in the beam formula or by relative phase
            % Here we return (n + m) * psi(z) to be used as an ADDITIVE shift to the Gaussian Gouy phase
            phi = (obj.n + obj.m) .* obj.GouyPhase;
        end
        
        function wH = get.HermiteWaistX(obj)
            % Standard Hermite spot size in X: w(z) * sqrt(n + 1)
            % Note: some conventions use sqrt(2n + 1). We'll follow the old code's intent.
            wH = obj.Waist .* sqrt(obj.n + 1);
        end
        
        function wH = get.HermiteWaistY(obj)
            wH = obj.Waist .* sqrt(obj.m + 1);
        end
        
        function wH = get.HermiteWaist(obj)
            % Combined waist (RMSE)
            wH = obj.Waist .* sqrt(obj.n + obj.m + 1);
        end
    end
    
    methods (Static)
        function wH = getWaistOneDirection(z, w0, zr, n)
            % Static method for one-dimensional Hermite waist
            w = w0 * sqrt(1 + (z/zr).^2);
            wH = w * sqrt(n + 1);
        end
        
        function wH = getWaist(z, w0, zr, n, m)
            % Static method for two-dimensional Hermite waist
            w = w0 * sqrt(1 + (z/zr).^2);
            wH = w * sqrt(n + m + 1);
        end
    end
end