classdef FFTUtils
    % FFTUtils - FFT utilities for beam propagation
    % Compatible with GNU Octave and MATLAB
    
    properties
        normalize = true;
        shiftFlag = true;
    end
    
    methods
        function obj = FFTUtils(normalize, shiftFlag)
            if nargin >= 1
                obj.normalize = normalize;
            end
            if nargin >= 2
                obj.shiftFlag = shiftFlag;
            end
        end
        
        function G = fft2(obj, g)
            if obj.shiftFlag
                G = fftshift(fft2(ifftshift(g)));
            else
                G = fft2(g);
            end
            
            if obj.normalize
                G = G / numel(g);
            end
        end
        
        function g = ifft2(obj, G)
            if obj.shiftFlag
                g = fftshift(ifft2(ifftshift(G)));
            else
                g = ifft2(G);
            end
            
            if obj.normalize
                g = g * numel(g);
            end
        end
        
        function G = fftn(obj, g)
            if obj.shiftFlag
                G = fftshift(fftn(ifftshift(g)));
            else
                G = fftn(g);
            end
            
            if obj.normalize
                G = G / numel(g);
            end
        end
        
        function g = ifftn(obj, G)
            if obj.shiftFlag
                g = fftshift(ifftn(ifftshift(G)));
            else
                g = ifftn(G);
            end
            
            if obj.normalize
                g = g * numel(g);
            end
        end
        
        function gprop = propagate(obj, g, kx, ky, z, lambda)
            G = obj.fft2(g);
            H = FFTUtils.transferFunction(kx, ky, z, lambda);
            gprop = obj.ifft2(G .* H);
        end
    end
    
    methods (Static)
        function H = transferFunction(kx, ky, z, lambda)
            k = 2*pi / lambda;
            kz = sqrt(k^2 - (kx.^2 + ky.^2));
            kz(real(kz) < 0) = 0;
            H = exp(1i * kz * z);
        end
        
        function H = transferSimple(kx, ky, z, lambda)
            k = 2*pi / lambda;
            phase = k*z - (kx.^2 + ky.^2)*z/(2*k);
            H = exp(1i * phase);
        end
        
        function G = fft2_centered(g)
            G = fftshift(fft2(ifftshift(g)));
        end
        
        function g = ifft2_centered(G)
            g = fftshift(ifft2(ifftshift(G)));
        end
        
        function G = fft2_std(g)
            G = fft2(g);
        end
        
        function g = ifft2_std(G)
            g = ifft2(G);
        end
    end
end
