classdef FFTUtils
%FFTUTILS Utility class for Fast Fourier Transform operations
%   This class provides normalized and consistent FFT operations for
%   beam propagation simulations using the angular spectrum method.
%
%   Example:
%       fftOps = FFTUtils;
%       G = fftOps.fft2(g);  % Forward FFT
%       g  = fftOps.ifft2(G); % Inverse FFT
    
    properties
        normalize          % Normalization flag (default: true)
        shiftFlag          % FFT shift flag (default: true)
    end
    
    methods
        function self = FFTUtils(normalize, shiftFlag)
            %FFTUTILS Constructor
            %   fftOps = FFTUtils() creates with defaults
            %   fftOps = FFTUtils(normalize, shiftFlag) customizes
            
            if nargin < 1
                self.normalize = true;
            else
                self.normalize = normalize;
            end
            
            if nargin < 2
                self.shiftFlag = true;
            else
                self.shiftFlag = shiftFlag;
            end
        end
        
        function G = fft2(self, g)
            %FFT2 2D Fast Fourier Transform
            %   G = fft2(g) performs forward FFT of field g
            %
            %   Uses centered FFT convention (fftshift before and after)
            %   for consistent zero-frequency at center
            
            if self.shiftFlag
                G = fftshift(fft2(ifftshift(g)));
            else
                G = fft2(g);
            end
            
            if self.normalize
                G = G / numel(g);
            end
        end
        
        function g = ifft2(self, G)
            %IFFT2 2D Inverse Fast Fourier Transform
            %   g = ifft2(G) performs inverse FFT of frequency field G
            
            if self.shiftFlag
                g = fftshift(ifft2(ifftshift(G)));
            else
                g = ifft2(G);
            end
            
            if self.normalize
                g = g * numel(g);
            end
        end
        
        function G = fftn(self, g)
            %FFTN N-Dimensional FFT
            %   G = fftn(g) performs forward FFT along all dimensions
            
            if self.shiftFlag
                G = fftshift(fftn(ifftshift(g)));
            else
                G = fftn(g);
            end
            
            if self.normalize
                G = G / numel(g);
            end
        end
        
        function g = ifftn(self, G)
            %IFFTN N-Dimensional Inverse FFT
            %   g = ifftn(G) performs inverse FFT along all dimensions
            
            if self.shiftFlag
                g = fftshift(ifftn(ifftshift(G)));
            else
                g = ifftn(G);
            end
            
            if self.normalize
                g = g * numel(g);
            end
        end
        
        function prop = transferFunction(self, kx, ky, z, wavelength, varargin)
            %TRANSFERFUNCTION Angular spectrum transfer function
            %   prop = transferFunction(kx, ky, z, lambda) returns H
            %   prop = transferFunction(kx, ky, z, lambda, n) includes
            %          refractive index
            %
            %   H = exp(i*k*z*sqrt(1 - (kx^2 + ky^2)/k^2))
            
            k = 2*pi / wavelength;
            
            if nargin >= 5
                n = varargin{1};
                k = k * n;
            else
                n = 1;
            end
            
            kz = sqrt(k^2 - (kx.^2 + ky.^2));
            kz(real(kz) < 0) = 0;  % Evanescent waves
            
            prop = exp(1i * kz * z);
        end
        
        function gprop = propagate(self, g, kx, ky, z, wavelength)
            %PROPAGATE Propagate field using angular spectrum method
            %   gprop = propagate(g, kx, ky, z, lambda) returns propagated
            %   field at distance z
            %
            %   Algorithm:
            %   1. FFT g -> G
            %   2. Multiply by transfer function H
            %   3. IFFT G*H -> gprop
            
            G = self.fft2(g);
            H = self.transferFunction(kx, ky, z, wavelength);
            gprop = self.ifft2(G .* H);
        end
        
        function setNormalize(self, flag)
            %SETNORMALIZE Toggle normalization
            self.normalize = flag;
        end
        
        function setShiftFlag(self, flag)
            %SETSHIFTFLAG Toggle FFT shift
            self.shiftFlag = flag;
        end
    end
    
    methods (Static)
        function g = fft2_std(g)
            %FFT2_STD Standard MATLAB fft2 (no shift)
            g = fft2(g);
        end
        
        function g = ifft2_std(g)
            %IFFT2_STD Standard MATLAB ifft2 (no shift)
            g = ifft2(g);
        end
        
        function G = fft2_centered(g)
            %FFT2_CENTERED Centered FFT (shifted)
            G = fftshift(fft2(ifftshift(g)));
        end
        
        function g = ifft2_centered(g)
            %IFFT2_CENTERED Centered inverse FFT
            g = fftshift(ifft2(ifftshift(g)));
        end
        
        function [kx, ky] = freqVectors(N, D)
            %FREQVECTORS Create frequency vectors from grid parameters
            %   [kx, ky] = freqVectors(N, D) returns wave numbers
            %
            %   Example:
            %       [kx, ky] = freqVectors(1024, 1e-3);
            
            n = -N/2:N/2-1;
            du = 1 / D;
            u = n * du;
            kx = 2*pi*u;
            ky = kx;
        end
        
        function H = transferSimple(kx, ky, z, lambda)
            %TRANSFERSIMPLE Simple transfer function (paraxial)
            %   H = transferSimple(kx, ky, z, lambda) paraxial approximation
            %
            %   H = exp(i*k*z*(1 - (kx^2 + ky^2)/(2*k^2)))
            
            k = 2*pi / lambda;
            phase = k*z - (kx.^2 + ky.^2)*z/(2*k);
            H = exp(1i * phase);
        end
    end
end
