classdef FFTUtils
%FFTUTILS Utility class for FFT operations
    
    properties
        normalize = true;
        shiftFlag = true;
    end
    
    methods
        function self = FFTUtils(normalize, shiftFlag)
            if nargin >= 1, self.normalize = normalize; end
            if nargin >= 2, self.shiftFlag = shiftFlag; end
        end
        
        function G = fft2(self, g)
            if self.shiftFlag
                G = fftshift(fft2(ifftshift(g)));
            else
                G = fft2(g);
            end
            if self.normalize, G = G / numel(g); end
        end
        
        function g = ifft2(self, G)
            if self.shiftFlag
                g = fftshift(ifft2(ifftshift(G)));
            else
                g = ifft2(G);
            end
            if self.normalize, g = g * numel(g); end
        end
        
        function prop = transferFunction(self, kx, ky, z, wavelength)
            k = 2*pi / wavelength;
            kz = sqrt(k^2 - (kx.^2 + ky.^2));
            kz(real(kz) < 0) = 0;
            prop = exp(1i * kz * z);
        end
        
        function gprop = propagate(self, g, kx, ky, z, wavelength)
            G = self.fft2(g);
            H = self.transferFunction(kx, ky, z, wavelength);
            gprop = self.ifft2(G .* H);
        end
    end
    
    methods (Static)
        function G = fft2_centered(g)
            G = fftshift(fft2(ifftshift(g)));
        end
        
        function g = ifft2_centered(g)
            g = fftshift(ifft2(ifftshift(g)));
        end
        
        function H = transferSimple(kx, ky, z, lambda)
            k = 2*pi / lambda;
            phase = k*z - (kx.^2 + ky.^2)*z/(2*k);
            H = exp(1i * phase);
        end
    end
end
