classdef Wavefront
    % Wavefront - Optical wavefront analysis class
    % Compatible with GNU Octave and MATLAB
    %
    % Provides wavefront extraction, Zernike decomposition, and metrics
    % for optical beam characterization.
    %
    % Usage:
    %   wf = Wavefront(E, lambda)              % minimal constructor
    %   wf = Wavefront(E, lambda, grid)       % with GridUtils
    %
    %   phi = wf.getPhase();                  % wrapped phase
    %   I = wf.getIntensity();                % |E|^2
    %   coeffs = wf.fitZernike(36);          % fit 36 Zernike terms
    %   rms = wf.computeRMS();               % wavefront RMS
    %   strehl = wf.computeStrehl();         % Strehl ratio

    properties
        Field               % Complex field [Ny x Nx]
        Lambda = 632.8e-9   % Wavelength [m]
        dx = 1e-5           % Grid spacing x [m]
        dy = 1e-5            % Grid spacing y [m]
        Dx = 1e-3            % Grid extent x [m]
        Dy = 1e-3            % Grid extent y [m]
    end

    properties (Dependent, Hidden)
        Ny, Nx              % Grid dimensions (derived from Field)
    end

    methods
        function obj = Wavefront(E, lambda, varargin)
            % Constructor
            %   wf = Wavefront(E, lambda)              % minimal
            %   wf = Wavefront(E, lambda, grid)         % with GridUtils

            if nargin == 0
                return; % empty object
            end

            obj.Field = E;
            obj.Lambda = lambda;

            % Handle optional grid argument
            if nargin >= 3
                grid = varargin{1};
                if isa(grid, 'GridUtils')
                    obj.dx = grid.dx;
                    obj.dy = grid.dy;
                    obj.Dx = grid.Dx;
                    obj.Dy = grid.Dy;
                else
                    error('Wavefront:invalidGrid', ...
                        'Third argument must be a GridUtils instance');
                end
            end
        end

        % -----------------------------------------------------------------
        % Dependent properties
        % -----------------------------------------------------------------
        function Ny = get.Ny(obj)
            Ny = size(obj.Field, 1);
        end

        function Nx = get.Nx(obj)
            Nx = size(obj.Field, 2);
        end

        % -----------------------------------------------------------------
        % Field getters
        % -----------------------------------------------------------------

        function E = getField(obj)
            % getField - Return the complex field
            E = obj.Field;
        end

        function I = getIntensity(obj)
            % getIntensity - Return intensity |E|^2
            I = abs(obj.Field).^2;
        end

        function phi = getPhase(obj)
            % getPhase - Return wrapped phase angle(E) in [-pi, pi]
            phi = angle(obj.Field);
        end

        % -----------------------------------------------------------------
        % Grid helpers
        % -----------------------------------------------------------------

        function [rho, theta] = gridPolar(obj)
            % gridPolar - Return polar grid coordinates (normalized)
            %   rho normalized to [0, 1] at edge of grid
            %   theta in [-pi, pi]

            [X, Y] = obj.gridCartesian();

            % Normalize to unit circle (use max extent)
            maxR = min(obj.Dx/2, obj.Dy/2);
            rho = sqrt(X.^2 + Y.^2) / maxR;
            theta = atan2(Y, X);

            % Clip rho to [0, 1]
            rho = min(1, max(0, rho));
        end

        function [X, Y] = gridCartesian(obj)
            % gridCartesian - Return Cartesian grid arrays
            nx = (-obj.Nx/2:obj.Nx/2-1) * obj.dx;
            ny = (-obj.Ny/2:obj.Ny/2-1) * obj.dy;
            [X, Y] = meshgrid(nx, ny);
        end

        % -----------------------------------------------------------------
        % Zernike fitting
        % -----------------------------------------------------------------

        function coeffs = fitZernike(obj, nTerms)
            % fitZernike - Fit Zernike polynomial coefficients (least-squares)
            %
            %   coeffs = wf.fitZernike(36)  % fit 36 terms
            %
            % Returns:
            %   coeffs - [nTerms x 1] column vector of Zernike coefficients
            %
            % Notes:
            %   - Phase is unwrapped before fitting to handle discontinuities
            %     from wrapped [-pi, pi] phase for large wavefront excursions.
            %   - Only in-pupil pixels (rho <= 1) are included in the fit;
            %     out-of-pupil pixels are masked out to avoid biasing the
            %     recovered coefficients.

            if nTerms < 1 || nTerms > 36
                error('Wavefront:invalidNTerms', ...
                    'nTerms must be 1-36, got %d', nTerms);
            end

            [rho, theta] = obj.gridPolar();

            % Build pupil mask (1 inside aperture, 0 outside)
            pupilMask = rho <= 1;
            inPupil = pupilMask(:);

            % Get wrapped phase, then unwrap for reliable Zernike decomposition
            phi_wrapped = obj.getPhase();
            % Unwrap along rows, then along columns to handle 2D phase maps
            phi_unwrapped = unwrap(unwrap(phi_wrapped, [], 2), [], 1);

            % Apply pupil mask and flatten
            phi_masked = phi_unwrapped(:);
            rho_masked = rho(:);
            theta_masked = theta(:);

            % Keep only in-pupil points
            phi_masked = phi_masked(inPupil);
            rho_masked = rho_masked(inPupil);
            theta_masked = theta_masked(inPupil);

            % Build design matrix and fit via pseudo-inverse
            M = ZernikeUtils.zernikeMatrix(rho_masked, theta_masked, nTerms);
            coeffs = pinv(M) * phi_masked;  % least-squares solution
            coeffs = coeffs(:);  % ensure column vector
        end

        function phi_recon = reconstructZernike(obj, coeffs)
            % reconstructZernike - Reconstruct wavefront from Zernike coefficients
            %
            %   phi_recon = wf.reconstructZernike(coeffs);

            [rho, theta] = obj.gridPolar();
            nTerms = length(coeffs);

            M = ZernikeUtils.zernikeMatrix(rho, theta, nTerms);
            phi_recon = reshape(M * coeffs, obj.Ny, obj.Nx);
        end

        function residual = zernikeResidual(obj, nTerms)
            % zernikeResidual - RMS residual after fitting nTerms Zernikes
            % Computed only over in-pupil pixels to avoid edge artifacts.
            coeffs = obj.fitZernike(nTerms);
            [rho, theta] = obj.gridPolar();
            pupilMask = rho <= 1;
            inPupil = pupilMask(:);

            phi_fit = obj.reconstructZernike(coeffs);
            phi_unwrapped = unwrap(unwrap(obj.getPhase(), [], 2), [], 1);

            phi_fit_masked = phi_fit(:);
            phi_orig_masked = phi_unwrapped(:);

            phi_fit_masked = phi_fit_masked(inPupil);
            phi_orig_masked = phi_orig_masked(inPupil);

            residual = sqrt(mean((phi_orig_masked - phi_fit_masked).^2));
        end

        % -----------------------------------------------------------------
        % Metrics
        % -----------------------------------------------------------------

        function rms = computeRMS(obj, varargin)
            % computeRMS - Root-mean-square wavefront error
            %
            %   rms = wf.computeRMS()              % RMS of full phase
            %   rms = wf.computeRMS(coeffs)        % RMS from fitted Zernikes only

            if nargin == 1
                % RMS of actual phase
                phi = obj.getPhase();
                rms = sqrt(mean(phi(:).^2));
            else
                % RMS from fitted coefficients only (excludes residual)
                coeffs = varargin{1};
                phi_fit = obj.reconstructZernike(coeffs);
                rms = sqrt(mean(phi_fit(:).^2));
            end
        end

        function pv = computePV(obj)
            % computePV - Peak-to-valley wavefront error
            phi = obj.getPhase();
            pv = max(phi(:)) - min(phi(:));
        end

        function strehl = computeStrehl(obj)
            % computeStrehl - Strehl ratio via Maréchal approximation
            %
            %   strehl = exp(-sigma^2)
            %
            % where sigma is the RMS wavefront error in radians (phase units).
            % This is the standard Maréchal approximation: for small aberrations,
            % the Strehl ratio is approximately exp(-(RMS wavefront error in waves)^2).
            % Since computeRMS returns phase in radians, squaring gives (sigma_rad)^2
            % which equals (sigma_waves * 2*pi)^2 = (2*pi * sigma_waves)^2, but the
            % exponent exp(-sigma_rad^2) is the correct form for phase RMS in radians.

            sigma = obj.computeRMS();
            strehl = exp(-sigma^2);

            % Clamp to [0, 1] (numerical precision can exceed bounds)
            strehl = max(0, min(1, strehl));
        end

        function metrics = getMetrics(obj, nTerms)
            % getMetrics - Compute all wavefront metrics
            %
            %   metrics = wf.getMetrics(nTerms)   % include Zernike fitting
            %   metrics = wf.getMetrics()          % use 36 terms

            if nargin < 2
                nTerms = 36;
            end

            coeffs = obj.fitZernike(nTerms);

            metrics = struct();
            metrics.rms = obj.computeRMS(coeffs);
            metrics.pv = obj.computePV();
            metrics.strehl = obj.computeStrehl();
            metrics.residualRMS = obj.zernikeResidual(nTerms);
            metrics.nTerms = nTerms;
            metrics.coeffs = coeffs;
        end

        % -----------------------------------------------------------------
        % Visualization
        % -----------------------------------------------------------------

        function plotWavefront(obj, varargin)
            % plotWavefront - Display 2D wavefront phase map
            %
            %   wf.plotWavefront()
            %   wf.plotWavefront('units', 'waves')  % show in wavelength units

            parseargs = inputParser();
            addParameter(parseargs, 'units', 'rad');
            parse(parseargs, varargin{:});

            phi = obj.getPhase();

            figure;
            imagesc(phi);
            colorbar;
            axis square;
            colormap('hsv');
            caxis([-pi, pi]);

            if strcmp(parseargs.Results.units, 'waves')
                phi_waves = phi / (2*pi);
                imagesc(phi_waves);
                colorbar;
                caxis([-0.5, 0.5]);
                ylabel('Phase [waves]');
            else
                ylabel('Phase [rad]');
            end

            xlabel('x'); ylabel('y');
            title('Wavefront Phase');
        end

        function plotIntensity(obj)
            % plotIntensity - Display 2D intensity map

            I = obj.getIntensity();

            figure;
            imagesc(I);
            colorbar;
            axis square;
            colormap('hot');
            xlabel('x'); ylabel('y');
            title('Intensity |E|^2');
        end

        function plotZernikeCoeffs(obj, coeffs)
            % plotZernikeCoeffs - Bar chart of Zernike coefficients
            %
            %   wf.plotZernikeCoeffs(coeffs);

            n = length(coeffs);
            names = arrayfun(@(i) ZernikeUtils.zernikeName(i), 1:n, 'UniformOutput', false);

            figure;
            bar(1:n, coeffs);
            set(gca, 'XTick', 1:n);
            set(gca, 'XTickLabel', names, 'XTickLabelRotation', 60);
            ylabel('Coefficient [rad]');
            title('Zernike Decomposition');
            grid on;
        end

        function plotPhaseSlice(obj, plane, idx)
            % plotPhaseSlice - 1D phase cross-section
            %
            %   wf.plotPhaseSlice('x', Ny/2)   % horizontal slice at center
            %   wf.plotPhaseSlice('y', Nx/2)   % vertical slice at center

            phi = obj.getPhase();

            if strcmpi(plane, 'x')
                slice = phi(idx, :);
                x = (-obj.Nx/2:obj.Nx/2-1) * obj.dx;
                xlabel_str = 'x [m]';
            else
                slice = phi(:, idx)';
                x = (-obj.Ny/2:obj.Ny/2-1) * obj.dy;
                xlabel_str = 'y [m]';
            end

            figure;
            plot(x, slice);
            xlabel(xlabel_str);
            ylabel('Phase [rad]');
            title(sprintf('Phase Slice (%s=%d)', plane, idx));
            grid on;
        end
    end
end