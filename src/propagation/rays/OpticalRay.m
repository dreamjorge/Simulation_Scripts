classdef OpticalRay
    % OpticalRay - Ray in Cartesian coordinates (x, y, z)
    % Used in Hermite-based ray-tracing propagation.
    %
    % Parallel class hierarchy to CylindricalRay which uses cylindrical
    % coordinates (r, theta, z).
    % Properties follow the naming convention expected by assignCoordinates2CartesianRay
    % and copyRay2ArrayRay/copyArrayRay2Ray functions.

    properties
        xCoordinate     = 0   % Cartesian x coordinate [m]
        yCoordinate     = 0   % Cartesian y coordinate [m]
        zCoordinate     = 0   % Axial coordinate [m]

        zxSlope         = Inf % Slope dz/dx (Inf at z=0: ray propagates axially)
        zySlope         = Inf % Slope dz/dy
        xySlope         = Inf % Slope dy/dx

        hankelType      = 1   % 1 or 2: selects H^(1) or H^(2) Hankel function
    end

    methods
        function obj = OpticalRay(varargin)
            % Default constructor creates an empty ray at the origin.
            % No arguments required; all properties have defaults.
        end
    end
end
