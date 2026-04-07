classdef GridUtils
    % GridUtils - Grid generation utilities
    % Compatible with GNU Octave and MATLAB
    
    properties
        Nx, Ny, Nz
        Dx, Dy, Dz
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
                    obj.dz = Dz / Nz;
                end
                obj.dx = Dx / Nx;
                obj.dy = Dy / Ny;
            end
        end
        
        function [X, Y] = create2DGrid(obj)
            nx = -obj.Nx/2:obj.Nx/2-1;
            ny = -obj.Ny/2:obj.Ny/2-1;
            x = nx * (obj.Dx / obj.Nx);
            y = ny * (obj.Dy / obj.Ny);
            [X, Y] = meshgrid(x, y);
        end
        
        function [Kx, Ky] = createFreqGrid(obj)
            % Ensure properties are updated before grid generation
            % as they might have been changed manually in the script
            obj.dx = obj.Dx / obj.Nx;
            obj.dy = obj.Dy / obj.Ny;
            
            nx = -obj.Nx/2:obj.Nx/2-1;
            ny = -obj.Ny/2:obj.Ny/2-1;
            dux = 1 / obj.Dx;
            duy = 1 / obj.Dy;
            u = nx * dux;
            v = ny * duy;
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
            [Kx, Ky] = meshgrid(2*pi*u);
        end
        
        function [r, theta] = polarGrid(N, D)
            [X, Y] = GridUtils.meshgrid2D(N, D);
            [r, theta] = cart2pol(X, Y);
        end
    end
end
