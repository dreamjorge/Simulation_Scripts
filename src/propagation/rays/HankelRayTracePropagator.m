classdef HankelRayTracePropagator < IPropagator
    % HankelRayTracePropagator - Hankel-aware ray tracing (Strategy pattern)
    % Compatible with GNU Octave and MATLAB
    %
    % Constructor:
    %   prop = HankelRayTracePropagator(grid, hankelType)
    %   prop = HankelRayTracePropagator(grid, hankelType, method, dz, z0)
    %
    % Usage:
    %   beam = HankelLaguerre(w0, lambda, l, p, 2);
    %   prop = HankelRayTracePropagator(grid, 2, 'RK4', 1e-4);
    %   bundle = prop.propagate(beam, z_final);
    %
    % Wraps HankelRayTracer.propagate() with the IPropagator interface.
    % Tracks per-ray Hankel type and handles axis-crossing for Laguerre beams.
    %
    % For beams without axis-crossing concerns (HankelHermite), this propagator
    % still works correctly — it simply maintains the initial Hankel type.

    properties
        Grid        % GridUtils object
        HankelType  % Initial Hankel type for all rays
        Method      % 'RK4' or 'Euler'
        dz          % Step size (m); [] = auto
        z0          % Initial z plane
    end

    methods
        function obj = HankelRayTracePropagator(grid, hankelType, method, dz, z0)
            obj.Grid = grid;
            if nargin >= 2 && ~isempty(hankelType)
                obj.HankelType = hankelType;
            else
                obj.HankelType = 1;
            end
            if nargin >= 3 && ~isempty(method)
                obj.Method = method;
            else
                obj.Method = 'RK4';
            end
            if nargin >= 4 && ~isempty(dz)
                obj.dz = dz;
            else
                obj.dz = [];
            end
            if nargin >= 5
                obj.z0 = z0;
            else
                obj.z0 = 0;
            end
        end

        function bundle = propagate(obj, beam, z_final)
            [X, Y] = obj.Grid.create2DGrid();

            if isempty(obj.dz)
                dz_step = (z_final - obj.z0) / 100;
            else
                dz_step = obj.dz;
            end

            bundle = RayBundle(X, Y, obj.z0);
            bundle.ht(:) = obj.HankelType;
            bundle = HankelRayTracer.propagate(bundle, beam, z_final, dz_step, obj.Method);
        end
    end
end
