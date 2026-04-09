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
    
    methods
        function field = opticalField(obj, X, Y, z)
            % This is a dummy implementation for portability.
            % Must be overridden by subclasses.
            error('ParaxialBeam:opticalFieldNotImplemented', ...
                'The method opticalField must be implemented by the subclass.');
        end
    end
end
