classdef IPropagator < handle
    % IPropagator - Abstract base class for beam propagators (Strategy pattern)
    % Compatible with GNU Octave and MATLAB
    %
    % All propagators implement a common interface:
    %
    %   result = prop.propagate(beam, z_final)
    %
    % where:
    %   beam    — any ParaxialBeam subclass
    %   z_final — target axial position (m)
    %   result  — propagator-dependent output (complex field matrix or RayBundle)
    %
    % This allows swapping propagation methods without modifying the beam or
    % the calling code — the canonical Strategy pattern.
    %
    % Concrete subclasses:
    %   FFTPropagator       — Angular spectrum method via FFT
    %   AnalyticPropagator  — Direct evaluation of beam.opticalField at z_final
    %   RayTracePropagator  — Phase-gradient ray tracing via RayTracer
    %
    % Usage example:
    %   beam = GaussianBeam(100e-6, 632.8e-9);
    %   grid = GridUtils(256, 256, 1e-3, 1e-3);
    %
    %   prop = FFTPropagator(grid, 632.8e-9);
    %   field = prop.propagate(beam, 0.1);
    %
    %   prop2 = RayTracePropagator(grid, 'RK4', 1e-3);
    %   bundle = prop2.propagate(beam, 0.1);

    methods
        function result = propagate(obj, beam, z_final)
            % propagate - Propagate beam to z_final.
            %
            % Subclasses MUST override this method.
            %
            % Parameters:
            %   beam    (ParaxialBeam): beam to propagate
            %   z_final (scalar):       target axial position (m)
            %
            % Returns:
            %   result: propagator-specific output.
            %           FFTPropagator       -> complex field matrix [Ny x Nx]
            %           AnalyticPropagator  -> complex field matrix [Ny x Nx]
            %           RayTracePropagator  -> RayBundle object
            error('IPropagator:notImplemented', ...
                '%s must implement propagate(obj, beam, z_final).', class(obj));
        end
    end
end
