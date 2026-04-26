classdef AnalyticPropagator < IPropagator
    % AnalyticPropagator - Direct analytical propagation via beam.opticalField
    % Compatible with GNU Octave and MATLAB
    %
    % Constructor:
    %   prop = AnalyticPropagator(grid)
    %   prop = AnalyticPropagator(grid, z0)
    %
    % Usage:
    %   field = prop.propagate(beam, z_final)
    %
    % Strategy:
    %   Directly evaluates beam.opticalField(X, Y, z_final) on the grid.
    %   This works for all ParaxialBeam subclasses that have a closed-form
    %   expression (Gaussian, Hermite-Gaussian, Laguerre-Gaussian, etc.).
    %
    %   No FFT or ray stepping is involved — each call is O(N^2).
    %
    % When to use:
    %   - Any time you trust the beam's analytic formula to be exact.
    %   - As the reference / ground truth for comparing FFT or ray results.
    %   - When diffraction effects beyond the paraxial approximation are negligible.

    properties
        Grid    % GridUtils object — defines the (X, Y) sample points
        z0      % Initial z position (default 0)
    end

    methods
        function obj = AnalyticPropagator(grid, z0)
            % Constructor
            % grid: GridUtils object
            % z0:   initial z plane (default 0)

            % Emit deprecation warning (Strangler Fig migration)
            warning('BeamFactory:deprecated', ...
                'src/propagation/field/AnalyticPropagator is deprecated. Use +paraxial/+propagation/+field/AnalyticPropagator directly.');

            obj.Grid = grid;
            if nargin >= 2
                obj.z0 = z0;
            else
                obj.z0 = 0;
            end
        end

        function field = propagate(obj, beam, z_final)
            % propagate - Evaluate beam field analytically at z_final.
            %
            % Parameters:
            %   beam    (ParaxialBeam): beam to evaluate
            %   z_final (scalar):       target z (m)
            %
            % Returns:
            %   field [Ny x Nx]: complex optical field at z_final
            [X, Y] = obj.Grid.create2DGrid();
            field  = beam.opticalField(X, Y, z_final);
        end
    end
end
