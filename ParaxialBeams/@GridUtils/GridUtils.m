classdef GridUtils
%GRIDUTILS Utility class for grid generation
    
    properties
        dx, dy, dz, Nx, Ny, Nz, Dx, Dy, Dz
    end
    
    methods
        function self = GridUtils(Nx, Ny, Dx, Dy, Nz, Dz)
            if nargin >= 4
                self.Nx = Nx; self.Ny = Ny;
                self.Dx = Dx; self.Dy = Dy;
                self.dx = Dx / Nx; self.dy = Dy / Ny;
            end
            if nargin >= 6
                self.Nz = Nz; self.Dz = Dz;
                self.dz = Dz / Nz;
            end
        end
        
        function [X, Y] = create2DGrid(self)
            n = -self.Nx/2:self.Nx/2-1;
            x = n * self.dx;
            [X, Y] = meshgrid(x, x);
        end
        
        function [Kx, Ky] = createFreqGrid(self)
            n = -self.Nx/2:self.Nx/2-1;
            du = 1 / self.Dx;
            u = n * du;
            [Kx, Ky] = meshgrid(2*pi*u);
        end
    end
    
    methods (Static)
        function [X, Y] = meshgrid2D(N, D)
            n = -N/2:N/2-1;
            x = n * (D / N);
            [X, Y] = meshgrid(x, x);
        end
        
        function [Kx, Ky] = freqGrid(N, D)
            n = -N/2:N/2-1;
            u = n * (1 / D);
            [Kx, Ky] = meshgrid(2*pi*u);
        end
    end
end
