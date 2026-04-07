classdef GridUtils
%GRIDUTILS Utility class for grid generation and management
%   This class provides methods for creating and managing computational
%   grids for beam propagation simulations.
%
%   Example:
%       grid = GridUtils;
%       [X, Y, Kx, Ky] = grid.create2DGrid(1024, 100e-6);
    
    properties (Hidden)
        dx          % Spatial resolution
        dy          % Spatial resolution (y-direction)
        dz          % Propagation resolution
        Nx          % Number of points x
        Ny          % Number of points y
        Nz          % Number of points z
        Dx          % Window size x
        Dy          % Window size y
        Dz          % Window size z
    end
    
    methods
        function self = GridUtils(Nx, Ny, Dx, Dy, Nz, Dz)
            %GRIDUTILS Constructor
            %   grid = GridUtils(Nx, Ny, Dx, Dy) creates a 2D grid
            %   grid = GridUtils(Nx, Ny, Dx, Dy, Nz, Dz) creates a 3D grid
            
            if nargin >= 4
                self.Nx = Nx;
                self.Ny = Ny;
                self.Dx = Dx;
                self.Dy = Dy;
                self.dx = Dx / Nx;
                self.dy = Dy / Ny;
            end
            
            if nargin >= 6
                self.Nz = Nz;
                self.Dz = Dz;
                self.dz = Dz / Nz;
            end
        end
        
        function [X, Y] = create2DGrid(self)
            %CREATE2DGRID Create 2D spatial grid
            %   [X, Y] = create2DGrid returns meshgrid arrays
            
            n = -self.Nx/2:self.Nx/2-1;
            x = n * self.dx;
            y = x;
            [X, Y] = meshgrid(x, y);
        end
        
        function [X, Y, Z] = create3DGrid(self)
            %CREATE3DGRID Create 3D spatial grid
            %   [X, Y, Z] = create3DGrid returns 3D meshgrid arrays
            
            n = -self.Nx/2:self.Nx/2-1;
            x = n * self.dx;
            y = x;
            z = (0:self.Nz-1) * self.dz;
            [X, Y, Z] = meshgrid(x, y, z);
        end
        
        function [Kx, Ky] = createFreqGrid(self)
            %CREATEFREQGRID Create frequency (k-space) grid
            %   [Kx, Ky] = createFreqGrid returns k-space vectors
            
            n = -self.Nx/2:self.Nx/2-1;
            du = 1 / self.Dx;
            u = n * du;
            [Kx, Ky] = meshgrid(2*pi*u);
        end
        
        function [U, V] = createNormFreqGrid(self)
            %CREATENORMFREQGRID Create normalized frequency grid
            %   [U, V] = createNormFreqGrid returns normalized u,v
            
            n = -self.Nx/2:self.Nx/2-1;
            du = 1 / self.Dx;
            U = n * du;
            V = U;
            [U, V] = meshgrid(U, V);
        end
        
        function [r, theta] = createPolarGrid(self)
            %CREATEPOLARGRID Create polar grid from Cartesian
            %   [r, theta] = createPolarGrid returns polar coordinates
            
            [X, Y] = self.create2DGrid();
            [r, theta] = cart2pol(X, Y);
        end
        
        function estimateResolution(self, maxWaist, safetyFactor)
            %ESTIMATERESOLUTION Estimate optimal resolution
            %   estimateResolution(maxWaist) calculates dx for given max waist
            %   estimateResolution(maxWaist, safetyFactor) uses custom safety
            
            if nargin < 2
                safetyFactor = 1.1;
            end
            
            self.Dx = safetyFactor * 2 * maxWaist;
            self.dx = self.Dx / self.Nx;
            self.Dy = self.Dx;
            self.dy = self.Dx / self.Nx;
        end
        
        function info = getGridInfo(self)
            %GETGRIDINFO Return grid parameters as struct
            
            info = struct(...
                'Nx', self.Nx, ...
                'Ny', self.Ny, ...
                'Dx', self.Dx, ...
                'Dy', self.Dy, ...
                'dx', self.dx, ...
                'dy', self.dy);
            
            if ~isempty(self.Nz)
                info.Nz = self.Nz;
                info.Dz = self.Dz;
                info.dz = self.dz;
            end
        end
    end
    
    methods (Static)
        function [X, Y] = meshgrid2D(N, D)
            %MESHGRID2D Quick 2D meshgrid creation
            %   [X, Y] = meshgrid2D(N, D) creates NxN grid with size D
            
            n = -N/2:N/2-1;
            dx = D / N;
            x = n * dx;
            [X, Y] = meshgrid(x, x);
        end
        
        function [Kx, Ky] = freqGrid(N, D)
            %FREQGRID Quick frequency grid
            %   [Kx, Ky] = freqGrid(N, D) returns k-space for N points, size D
            
            n = -N/2:N/2-1;
            du = 1 / D;
            u = n * du;
            [Kx, Ky] = meshgrid(2*pi*u);
        end
        
        function [r, theta] = polarGrid(N, D)
            %POLARGRID Quick polar grid
            %   [r, theta] = polarGrid(N, D) returns polar coords
            
            [X, Y] = GridUtils.meshgrid2D(N, D);
            [r, theta] = cart2pol(X, Y);
        end
    end
end
