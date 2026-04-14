classdef RayBundle < handle
    % RayBundle Class to manage a set of optical rays.
    %
    % Properties:
    %   x, y, z - Coordinates (m)
    %   r, theta - Polar coordinates (rad, m)
    %   slopes - Local slopes (dx/dz, dy/dz, etc.)
    
    properties
        x % x-coordinates (Ny x Nx x Nz)
        y % y-coordinates (Ny x Nx x Nz)
        z % z-coordinates (Ny x Nx x Nz)
        r % radial coordinates
        theta % angular coordinates
        sx % dx/dz slopes
        sy % dy/dz slopes
        sz % dz/dz (always 1 for standard paraxial, but tracking for generality)
    end
    
    properties (Dependent)
        Nx
        Ny
        Nz
    end
    
    methods
        function obj = RayBundle(x0, y0, z0)
            % Initialize a bundle of rays at a starting z plane
            % x0, y0 are matrices/arrays of initial coordinates
            if nargin > 0
                obj.x = x0;
                obj.y = y0;
                obj.z = z0 * ones(size(x0));
                [obj.theta, obj.r] = cart2pol(obj.x, obj.y);
                obj.sx = zeros(size(x0));
                obj.sy = zeros(size(x0));
                obj.sz = ones(size(x0));
            end
        end
        
        function val = get.Nx(obj)
            val = size(obj.x, 2);
        end
        function val = get.Ny(obj)
            val = size(obj.x, 1);
        end
        function val = get.Nz(obj)
            val = size(obj.x, 3);
        end
        
        function addStep(obj, x, y, z, sx, sy)
            % Add a new slice of coordinates to the bundle
            obj.x = cat(3, obj.x, x);
            obj.y = cat(3, obj.y, y);
            obj.z = cat(3, obj.z, z);
            obj.sx = cat(3, obj.sx, sx);
            obj.sy = cat(3, obj.sy, sy);
        end
    end
    
    methods (Static)
        function bundle = createGrid(Nx, Ny, Dx, Dy, z0)
            if nargin < 5, z0 = 0; end
            [X, Y] = meshgrid(linspace(-Dx/2, Dx/2, Nx), linspace(-Dy/2, Dy/2, Ny));
            bundle = RayBundle(X, Y, z0);
        end
        
        function bundle = createConcentric(Nr, Ntheta, maxR, z0)
            if nargin < 4, z0 = 0; end
            r_vals = linspace(0, maxR, Nr);
            th_vals = linspace(0, 2*pi, Ntheta+1);
            th_vals(end) = []; % avoid duplicate at 0/2pi
            [R, TH] = meshgrid(r_vals, th_vals);
            [X, Y] = pol2cart(TH, R);
            bundle = RayBundle(X, Y, z0);
        end
        
        function bundle = createCrosshairs(N, L, z0)
            if nargin < 3, z0 = 0; end
            vals = linspace(-L/2, L/2, N);
            x_line = [vals, zeros(1, N)];
            y_line = [zeros(1, N), vals];
            bundle = RayBundle(x_line, y_line, z0);
        end
    end
end
