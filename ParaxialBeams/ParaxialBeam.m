classdef (Abstract) ParaxialBeam < handle
    % ParaxialBeam - Abstract base class for all paraxial beam models.
    % Enforces a unified API for optical field computation.
    
    properties
        Lambda % Wavelength (m)
        k % Wave number (rad/m)
    end
    
    methods
        function obj = ParaxialBeam(lambda)
            if nargin > 0
                obj.Lambda = lambda;
                obj.k = 2 * pi / lambda;
            end
        end
    end
    
    methods (Abstract)
        % field = opticalField(obj, X, Y, z)
        % Computes the complex optical field at given coordinates.
        % X, Y: matrix of coordinates (m)
        % z: scalar propagation distance (m)
        field = opticalField(obj, X, Y, z)
    end
end
