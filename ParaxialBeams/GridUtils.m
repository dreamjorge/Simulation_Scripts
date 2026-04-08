classdef GridUtils
    % GridUtils - Grid generation utilities
    % Compatible with GNU Octave and MATLAB
    
    properties
        Nx, Ny, Nz = 0
        Dx, Dy, Dz = 0
    end

    properties (Dependent)
        dx, dy, dz
    end

    methods
        function obj = GridUtils(Nx, Ny, Dx, Dy, Nz, Dz)
            if nargin > 0
                obj.Nx = Nx;
                obj.Ny = Ny;
                obj.Dx = Dx;
                obj.Dy = Dy;
                if nargin >= 6
                    obj.Nz = Nz;
                    obj.Dz = Dz;
                end
            end
        end

        function val = get.dx(obj)
            val = obj.Dx / obj.Nx;
        end

        function val = get.dy(obj)
            val = obj.Dy / obj.Ny;
        end

        function val = get.dz(obj)
            if isempty(obj.Nz) || obj.Nz == 0
                val = 0;
            else
                val = obj.Dz / obj.Nz;
            end
        end

        function [X, Y] = create2DGrid(obj)
            nx = -obj.Nx/2:obj.Nx/2-1;
            ny = -obj.Ny/2:obj.Ny/2-1;
            x = nx * obj.dx;
            y = ny * obj.dy;
            [X, Y] = meshgrid(x, y);
        end

        function [Kx, Ky] = createFreqGrid(obj)
            nx = -obj.Nx/2:obj.Nx/2-1;
            ny = -obj.Ny/2:obj.Ny/2-1;
            u = nx * (1 / obj.Dx);
            v = ny * (1 / obj.Dy);
            [Kx, Ky] = meshgrid(2*pi*u, 2*pi*v);
        end
        
        function [X, Y, Z] = create3DGrid(obj)
            nx = -obj.Nx/2:obj.Nx/2-1;
            ny = -obj.Ny/2:obj.Ny/2-1;
            nz = 0:obj.Nz-1;
            x = nx * (obj.Dx / obj.Nx);
            y = ny * (obj.Dy / obj.Ny);
            z = nz * (obj.Dz / obj.Nz);
            [X, Y, Z] = meshgrid(x, y, z);
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
            [Kx, Ky] = meshgrid(2*pi*u, 2*pi*u);
        end
        
        function [r, theta] = polarGrid(N, D)
            [X, Y] = GridUtils.meshgrid2D(N, D);
            [r, theta] = cart2pol(X, Y);
        end
    end
end
