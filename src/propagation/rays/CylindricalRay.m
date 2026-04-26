classdef CylindricalRay
    % CylindricalRay - Ray in cylindrical coordinates (r, theta, z)
    % Used in Hankel-based Laguerre beam ray-tracing propagation.
    %
    % Replaces the legacy @CylindricalRay class folder (removed in refactor).
    % Properties follow the same naming convention as the original class.

    properties
        rCoordinate     = 0   % Radial coordinate [m]
        thetaCoordinate = 0   % Angular coordinate [rad]
        zCoordinate     = 0   % Axial coordinate [m]
        xCoordinate     = 0   % Cartesian x (derived from r, theta)
        yCoordinate     = 0   % Cartesian y (derived from r, theta)

        zrSlope         = Inf % Slope dz/dr (Inf at z=0: ray propagates axially)
        zthSlope        = Inf % Slope dz/dtheta
        rthSlope        = Inf % Slope dr/dtheta

        hankelType      = 1   % 1 or 2: selects H^(1) or H^(2) Hankel function
    end

    methods
        function obj = CylindricalRay(varargin)
            % Default constructor creates an empty ray at the origin.
            % No arguments required; all properties have defaults.

            % Emit deprecation warning (Strangler Fig migration)
            warning('BeamFactory:deprecated', ...
                'src/propagation/rays/CylindricalRay is deprecated. Use +paraxial/+propagation/+rays/CylindricalRay directly.');
        end
    end
end
