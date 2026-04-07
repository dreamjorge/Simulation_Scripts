% PhysicalConstants - Physical constants and utility functions
% Compatible with GNU Octave and MATLAB
%
% Usage:
%   zr = PhysicalConstants.rayleighDistance(w0, lambda)
%   k = PhysicalConstants.waveNumber(lambda)

function varargout = PhysicalConstants(action, varargin)
    switch action
        case 'speed_of_light'
            varargout{1} = 299792458;
        case 'planck'
            varargout{1} = 6.62607015e-34;
        case 'planck_reduced'
            varargout{1} = 1.054571817e-34;
        case 'vacuum_permittivity'
            varargout{1} = 8.8541878128e-12;
        case 'vacuum_permeability'
            varargout{1} = 1.25663706212e-6;
        case 'impedance_vacuum'
            varargout{1} = 376.730313668;
        case 'waveNumber'
            lambda = varargin{1};
            varargout{1} = 2*pi ./ lambda;
        case 'rayleighDistance'
            w0 = varargin{1};
            lambda = varargin{2};
            varargout{1} = pi * w0.^2 ./ lambda;
        case 'waistAtZ'
            w0 = varargin{1};
            z = varargin{2};
            lambda = varargin{3};
            if nargin >= 5
                zr = varargin{4};
            else
                zr = pi * w0.^2 ./ lambda;
            end
            varargout{1} = w0 .* sqrt(1 + (z ./ zr).^2);
        case 'radiusOfCurvature'
            z = varargin{1};
            zr = varargin{2};
            varargout{1} = z .* (1 + (zr ./ z).^2);
        case 'gouyPhase'
            z = varargin{1};
            zr = varargin{2};
            varargout{1} = atan(z ./ zr);
        otherwise
            error('Unknown action: %s', action);
    end
end

% Shorthand functions for convenience
function zr = rayleighDistance(w0, lambda)
    zr = pi * w0.^2 ./ lambda;
end

function k = waveNumber(lambda)
    k = 2*pi ./ lambda;
end

function w = waistAtZ(w0, z, lambda, varargin)
    if nargin >= 4
        zr = varargin{1};
    else
        zr = pi * w0.^2 ./ lambda;
    end
    w = w0 .* sqrt(1 + (z ./ zr).^2);
end

function R = radiusOfCurvature(z, zr)
    R = z .* (1 + (zr ./ z).^2);
end

function gouy = gouyPhase(z, zr)
    gouy = atan(z ./ zr);
end
