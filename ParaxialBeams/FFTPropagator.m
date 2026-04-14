classdef FFTPropagator < IPropagator
    % FFTPropagator - Angular spectrum propagation via FFT transfer function
    % Compatible with GNU Octave and MATLAB
    %
    % Constructor:
    %   prop = FFTPropagator(grid, lambda)
    %   prop = FFTPropagator(grid, lambda, z0)
    %
    % Usage:
    %   field = prop.propagate(beam, z_final)
    %
    % Strategy (Angular Spectrum Method):
    %   1. Evaluate the beam field at z0 (the initial plane): u0 = beam.opticalField(X,Y,z0)
    %   2. Compute the 2D FFT: U0 = fft2(u0)
    %   3. Apply the free-space transfer function:
    %        H(kx, ky) = exp(i * kz * dz),   kz = sqrt(k^2 - kx^2 - ky^2)
    %      where dz = z_final - z0
    %   4. Return to spatial domain: u = ifft2(U0 .* H)
    %
    % Spatial frequencies (kx, ky) are built from the GridUtils frequency grid.
    %
    % When to use:
    %   - Full diffraction simulation including interference and evanescent effects.
    %   - When the beam profile is arbitrary (measured, SLM-shaped, etc.).
    %   - Benchmark: compare with AnalyticPropagator for paraxial beams.

    properties
        Grid    % GridUtils object — defines spatial and frequency sampling
        Lambda  % Wavelength (m)
        z0      % Initial z plane (default 0)
        FFT     % FFTUtils instance (normalize=true, shiftFlag=true)
    end

    methods
        function obj = FFTPropagator(grid, lambda, z0)
            % Constructor
            % grid:   GridUtils object
            % lambda: wavelength (m)
            % z0:     initial z plane (default 0)
            obj.Grid   = grid;
            obj.Lambda = lambda;
            obj.FFT    = FFTUtils();
            if nargin >= 3
                obj.z0 = z0;
            else
                obj.z0 = 0;
            end
        end

        function field = propagate(obj, beam, z_final)
            % propagate - Propagate beam from z0 to z_final via angular spectrum.
            %
            % Parameters:
            %   beam    (ParaxialBeam): beam to propagate
            %   z_final (scalar):       target z (m)
            %
            % Returns:
            %   field [Ny x Nx]: complex optical field at z_final
            [X, Y]   = obj.Grid.create2DGrid();
            [Kx, Ky] = obj.Grid.createFreqGrid();

            dz     = z_final - obj.z0;
            u0     = beam.opticalField(X, Y, obj.z0);
            field  = obj.FFT.propagate(u0, Kx, Ky, dz, obj.Lambda);
        end
    end
end
