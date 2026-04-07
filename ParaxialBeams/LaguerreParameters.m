classdef LaguerreParameters < GaussianParameters
    % LaguerreParameters - Parameters for Laguerre-Gaussian beams
    % Inherits from GaussianParameters and adds Laguerre-specific properties
    
    properties
        l % angular number (topological charge)
        p % radial number
    end
    
    properties (Dependent)
        LaguerreWaist  % Spot size of Laguerre-Gaussian Beam
        PhiPhase       % Total phase shift
    end
    
    methods
        function obj = LaguerreParameters(z, w0, lambda, l, p)
            % Constructor
            % z: propagation distance
            % w0: initial waist
            % lambda: wavelength
            % l: angular mode number
            % p: radial mode number
            
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
        
        %% Dependent property getters
        function phi = get.PhiPhase(obj)
            % Total phase shift: (abs(l) + 2*p) * GouyPhase
            % This is the additional phase shift relative to the fundamental Gaussian mode
            phi = (abs(obj.l) + 2*obj.p) .* obj.GouyPhase;
        end
        
        function wL = get.LaguerreWaist(obj)
            % Spot size of LG beam: w(z) * sqrt(2*p + abs(l) + 1)
            wL = obj.Waist .* sqrt(2*obj.p + abs(obj.l) + 1);
        end
    end
    
    methods (Static)
        function wL = getWaist(z, w0, zr, l, p)
            % Static method for Laguerre waist
            w = w0 * sqrt(1 + (z/zr).^2);
            wL = w * sqrt(2*p + abs(l) + 1);
        end
        
        function L = getAssociatedLaguerrePolynomial(p, l, x)
            % Calculates the associated Laguerre polynomial L_p^|l|(x)
            % Using the explicit summation formula
            n = p;
            k = abs(l);
            L = zeros(size(x));
            for m = 0:n
                term = ((-1)^m * nchoosek(n + k, n - m) .* x.^m) ./ factorial(m);
                L = L + term;
            end
        end
    end
end