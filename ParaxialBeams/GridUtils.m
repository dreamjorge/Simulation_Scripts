% GridUtils - Grid generation utilities
% Compatible with GNU Octave and MATLAB
%
% Usage:
%   grid = GridUtils(Nx, Ny, Dx, Dy)
%   [X, Y] = GridUtils.create2DGrid(Nx, Ny, Dx, Dy)
%   [Kx, Ky] = GridUtils.createFreqGrid(Nx, Ny, Dx, Dy)

function varargout = GridUtils(action, varargin)
    switch action
        case 'create2DGrid'
            Nx = varargin{1};
            Ny = varargin{2};
            Dx = varargin{3};
            Dy = varargin{4};
            
            nx = -Nx/2:Nx/2-1;
            ny = -Ny/2:Ny/2-1;
            x = nx * (Dx / Nx);
            y = ny * (Dy / Ny);
            [X, Y] = meshgrid(x, y);
            varargout{1} = X;
            varargout{2} = Y;
            
        case 'createFreqGrid'
            Nx = varargin{1};
            Ny = varargin{2};
            Dx = varargin{3};
            Dy = varargin{4};
            
            nx = -Nx/2:Nx/2-1;
            ny = -Ny/2:Ny/2-1;
            dux = 1 / Dx;
            duy = 1 / Dy;
            u = nx * dux;
            v = ny * duy;
            [Kx, Ky] = meshgrid(2*pi*u, 2*pi*v);
            varargout{1} = Kx;
            varargout{2} = Ky;
            
        case 'create3DGrid'
            Nx = varargin{1};
            Ny = varargin{2};
            Nz = varargin{3};
            Dx = varargin{4};
            Dy = varargin{5};
            Dz = varargin{6};
            
            nx = -Nx/2:Nx/2-1;
            ny = -Ny/2:Ny/2-1;
            nz = 0:Nz-1;
            x = nx * (Dx / Nx);
            y = ny * (Dy / Ny);
            z = nz * (Dz / Nz);
            [X, Y, Z] = meshgrid(x, y, z);
            varargout{1} = X;
            varargout{2} = Y;
            varargout{3} = Z;
            
        case 'meshgrid2D'
            N = varargin{1};
            D = varargin{2};
            n = -N/2:N/2-1;
            x = n * (D / N);
            [X, Y] = meshgrid(x, x);
            varargout{1} = X;
            varargout{2} = Y;
            
        case 'freqGrid'
            N = varargin{1};
            D = varargin{2};
            n = -N/2:N/2-1;
            u = n * (1 / D);
            [Kx, Ky] = meshgrid(2*pi*u);
            varargout{1} = Kx;
            varargout{2} = Ky;
            
        case 'polarGrid'
            N = varargin{1};
            D = varargin{2};
            [X, Y] = GridUtils('meshgrid2D', N, D);
            [R, Theta] = cart2pol(X, Y);
            varargout{1} = R;
            varargout{2} = Theta;
            
        otherwise
            error('Unknown action: %s', action);
    end
end

% Shorthand functions
function [X, Y] = create2DGrid(Nx, Ny, Dx, Dy)
    [X, Y] = GridUtils('create2DGrid', Nx, Ny, Dx, Dy);
end

function [Kx, Ky] = createFreqGrid(Nx, Ny, Dx, Dy)
    [Kx, Ky] = GridUtils('createFreqGrid', Nx, Ny, Dx, Dy);
end

function [X, Y, Z] = create3DGrid(Nx, Ny, Nz, Dx, Dy, Dz)
    [X, Y, Z] = GridUtils('create3DGrid', Nx, Ny, Nz, Dx, Dy, Dz);
end
