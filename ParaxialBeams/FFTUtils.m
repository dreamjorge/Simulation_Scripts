% FFTUtils - FFT utilities for beam propagation
% Compatible with GNU Octave and MATLAB
%
% Usage:
%   fftOps = FFTUtils()
%   G = FFTUtils.fft2(g, normalize, shiftFlag)
%   g = FFTUtils.ifft2(G, normalize, shiftFlag)
%   H = FFTUtils.transferFunction(kx, ky, z, lambda)

function varargout = FFTUtils(action, varargin)
    switch action
        case 'fft2'
            g = varargin{1};
            normalize = 2; % default: true
            shiftFlag = 1;
            if nargin >= 3
                normalize = varargin{2};
            end
            if nargin >= 4
                shiftFlag = varargin{3};
            end
            
            if shiftFlag
                G = fftshift(fft2(ifftshift(g)));
            else
                G = fft2(g);
            end
            
            if normalize
                G = G / numel(g);
            end
            varargout{1} = G;
            
        case 'ifft2'
            G = varargin{1};
            normalize = 2; % default: true
            shiftFlag = 1;
            if nargin >= 3
                normalize = varargin{2};
            end
            if nargin >= 4
                shiftFlag = varargin{3};
            end
            
            if shiftFlag
                g = fftshift(ifft2(ifftshift(G)));
            else
                g = ifft2(G);
            end
            
            if normalize
                g = g * numel(g);
            end
            varargout{1} = g;
            
        case 'fftn'
            g = varargin{1};
            normalize = 1;
            shiftFlag = 1;
            if nargin >= 3
                normalize = varargin{2};
            end
            if nargin >= 4
                shiftFlag = varargin{3};
            end
            
            if shiftFlag
                G = fftshift(fftn(ifftshift(g)));
            else
                G = fftn(g);
            end
            
            if normalize
                G = G / numel(g);
            end
            varargout{1} = G;
            
        case 'ifftn'
            G = varargin{1};
            normalize = 1;
            shiftFlag = 1;
            if nargin >= 3
                normalize = varargin{2};
            end
            if nargin >= 4
                shiftFlag = varargin{3};
            end
            
            if shiftFlag
                g = fftshift(ifftn(ifftshift(G)));
            else
                g = ifftn(G);
            end
            
            if normalize
                g = g * numel(g);
            end
            varargout{1} = g;
            
        case 'transferFunction'
            kx = varargin{1};
            ky = varargin{2};
            z = varargin{3};
            lambda = varargin{4};
            
            k = 2*pi / lambda;
            kz = sqrt(k^2 - (kx.^2 + ky.^2));
            kz(real(kz) < 0) = 0;
            H = exp(1i * kz * z);
            varargout{1} = H;
            
        case 'transferSimple'
            kx = varargin{1};
            ky = varargin{2};
            z = varargin{3};
            lambda = varargin{4};
            
            k = 2*pi / lambda;
            phase = k*z - (kx.^2 + ky.^2)*z/(2*k);
            H = exp(1i * phase);
            varargout{1} = H;
            
        case 'propagate'
            g = varargin{1};
            kx = varargin{2};
            ky = varargin{3};
            z = varargin{4};
            lambda = varargin{5};
            
            [G, normalize, shiftFlag] = deal(0, 1, 1);
            if nargin >= 7
                normalize = varargin{6};
            end
            if nargin >= 8
                shiftFlag = varargin{7};
            end
            
            G = FFTUtils('fft2', g, normalize, shiftFlag);
            H = FFTUtils('transferFunction', kx, ky, z, lambda);
            gprop = FFTUtils('ifft2', G .* H, normalize, shiftFlag);
            varargout{1} = gprop;
            
        case 'fft2_centered'
            g = varargin{1};
            G = fftshift(fft2(ifftshift(g)));
            varargout{1} = G;
            
        case 'ifft2_centered'
            g = varargin{1};
            g = fftshift(ifft2(ifftshift(g)));
            varargout{1} = g;
            
        case 'fft2_std'
            g = varargin{1};
            varargout{1} = fft2(g);
            
        case 'ifft2_std'
            g = varargin{1};
            varargout{1} = ifft2(g);
            
        otherwise
            error('Unknown action: %s', action);
    end
end

% Shorthand functions (only use when called from different context)
% These trigger the main function when called without 'action'
% Note: In Octave, nested functions can call the parent by name

% No shorthand functions to avoid recursion
% Use: G = FFTUtils('fft2', g, 1, 1);
