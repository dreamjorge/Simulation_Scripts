% GaussianParameters - Gaussian beam parameters
% Compatible with GNU Octave and MATLAB
%
% Usage:
%   params = GaussianParameters(z, w0, lambda)
%   zr = GaussianParameters('rayleighDistance', z, w0, lambda)

function varargout = GaussianParameters(action, varargin)
    if nargin < 1
        error('Usage: params = GaussianParameters(z, w0, lambda)');
    end
    
    % Direct method calls
    switch action
        case 'rayleighDistance'
            w0 = varargin{2};
            lambda = varargin{3};
            varargout{1} = pi * w0^2 ./ lambda;
            return;
        case 'getWaist'
            z = varargin{1};
            w0 = varargin{2};
            zr = varargin{3};
            varargout{1} = w0 * sqrt(1 + (z/zr)^2);
            return;
        case 'getPhase'
            z = varargin{1};
            zr = varargin{2};
            varargout{1} = atan(z/zr);
            return;
        case 'getRadius'
            z = varargin{1};
            zr = varargin{2};
            varargout{1} = z .* (1 + (zr ./ z).^2);
            return;
    end
    
    % Constructor: create params struct
    z = action;
    w0 = varargin{1};
    lambda = varargin{2};
    
    % Calculate dependent values (element-wise for vector z)
    zr = pi * w0^2 ./ lambda;
    k = 2*pi ./ lambda;
    w = w0 .* sqrt(1 + (z./zr).^2);
    gouy = atan(z ./ zr);
    R = z .* (1 + (zr ./ z).^2);
    amp = 1 ./ w;
    theta = atan(w0 ./ zr);
    
    % Create struct with function handles
    params = struct();
    params.zCoordinate = z;
    params.InitialWaist = w0;
    params.Wavelength = lambda;
    params.RayleighDistance = zr;
    params.k = k;
    params.Waist = w;
    params.GouyPhase = gouy;
    params.Radius = R;
    params.Amplitude = amp;
    params.DivergenceAngle = theta;
    
    params.toString = @() sprintf(...
        'GaussianParameters:\n  zCoordinate: %g\n  InitialWaist: %g\n  Wavelength: %g\n  RayleighDistance: %g\n  k: %g\n  Waist: %g\n', ...
        z, w0, lambda, zr, k, w);
    
    params.isEqual = @(other) ...
        abs(z - other.zCoordinate) < 1e-12 && ...
        abs(w0 - other.InitialWaist) < 1e-12 && ...
        abs(lambda - other.Wavelength) < 1e-12;
    
    varargout{1} = params;
end
