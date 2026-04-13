classdef AnalysisUtils
    % AnalysisUtils - Physical analysis and ray-tracing utilities
    % Consolidates logic for gradients, ray slopes, and wave combinations.
    
    methods (Static)
        function mzr = gradientRZ(fr, fz, k, dx, dz, x, z)
            % gradientRZ - Calculate ray slope in r-z plane
            % fr: Field slice at fixed r
            % fz: Field slice at fixed z
            % k: Wave number
            % dx, dz: Step sizes
            % x, z: Coordinates of interest
            
            % partial derivatives
            gz = gradient(fr) / dx;
            gr = gradient(fz) / dz + k;
            
            N = size(gr, 2);
            % slopes (using floor for coordinate lookup as in legacy code)
            idxZ = floor(z / dz) + 1;
            idxR = floor(x / dx) + floor(N / 2) + 1;
            
            mzr = gz(idxZ) / gr(idxR);
        end
        
        function [mzx, mzy, mxy] = gradientXYZ(fyz, fxz, fxy, k, dx, dy, dz, x, y, z)
            % gradientXYZ - Calculate local gradients in 3D
            % Consolidates local slope information for ray tracing.
            
            gx = gradient(fyz) / dx;
            gy = gradient(fxz) / dy;
            gz = gradient(fxy) / dz + k;
            
            % Lookup indices (adapted from legacy logic)
            idxX = floor(x / dx) + 1;
            idxY = floor(y / dy) + 1;
            idxZ = floor(z / dz) + 1;
            
            mzx = gx(idxX) / gz(idxZ);
            mzy = gy(idxY) / gz(idxZ);
            mxy = gx(idxX) / gy(idxY);
        end
        
    end
end
