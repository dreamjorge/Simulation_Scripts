classdef ZernikeUtils
    % ZernikeUtils - Static Zernike polynomial utilities (Noll indexing)
    % Compatible with GNU Octave and MATLAB
    %
    % Usage:
    %   Z = ZernikeUtils.zernike(n, rho, theta)   % compute Z_n
    %   name = ZernikeUtils.zernikeName(n)         % get name string
    %   M = ZernikeUtils.zernikeMatrix(rho, theta, N) % build design matrix
    %
    % Reference: Noll, R. J. (1976). Zernike polynomials and atmospheric turbulence.
    %   J. Opt. Soc. Am., 66(3), 207-211.

    methods (Static)
        function Z = zernike(n, rho, theta)
            % zernike - Compute Noll Zernike polynomial Z_n(rho, theta)
            %
            % Parameters:
            %   n     - Zernike index (Noll order, 1-based)
            %   rho   - radial coordinate [0, 1] (normalized)
            %   theta - azimuthal angle [rad]
            %
            % Returns:
            %   Z     - Zernike polynomial value(s)

            % Clamp rho to [0, 1] range
            rho = max(0, min(1, rho));

            switch n
                case 1
                    % Z1: Piston (constant)
                    Z = ones(size(rho));
                case 2
                    % Z2: Tilt X
                    Z = 2 * rho .* cos(theta);
                case 3
                    % Z3: Tilt Y
                    Z = 2 * rho .* sin(theta);
                case 4
                    % Z4: Defocus
                    Z = sqrt(3) .* (2 * rho.^2 - 1);
                case 5
                    % Z5: Astigmatism 0 degrees
                    Z = sqrt(6) * rho.^2 .* cos(2 * theta);
                case 6
                    % Z6: Astigmatism 45 degrees
                    Z = sqrt(6) * rho.^2 .* sin(2 * theta);
                case 7
                    % Z7: Coma X
                    Z = sqrt(8) * (3 * rho.^2 - 2) .* rho .* cos(theta);
                case 8
                    % Z8: Coma Y
                    Z = sqrt(8) * (3 * rho.^2 - 2) .* rho .* sin(theta);
                case 9
                    % Z9: Spherical aberration
                    Z = sqrt(8) * (6 * rho.^4 - 6 * rho.^2 + 1);
                case 10
                    % Z10: Secondary astigmatism 0 deg
                    Z = sqrt(10) * (4 * rho.^4 - 3 * rho.^2) .* cos(2 * theta);
                case 11
                    % Z11: Secondary astigmatism 45 deg
                    Z = sqrt(10) * (4 * rho.^4 - 3 * rho.^2) .* sin(2 * theta);
                case 12
                    % Z12: Secondary coma X
                    Z = sqrt(12) * (10 * rho.^6 - 12 * rho.^4 + 3 * rho.^2) .* rho .* cos(theta);
                case 13
                    % Z13: Secondary coma Y
                    Z = sqrt(12) * (10 * rho.^6 - 12 * rho.^4 + 3 * rho.^2) .* rho .* sin(theta);
                case 14
                    % Z14: Secondary spherical
                    Z = sqrt(12) * (20 * rho.^6 - 30 * rho.^4 + 12 * rho.^2 - 1);
                case 15
                    % Z15: Trefoil X
                    Z = sqrt(12) * (5 * rho.^6 - 4 * rho.^4) .* cos(3 * theta);
                case 16
                    % Z16: Trefoil Y
                    Z = sqrt(12) * (5 * rho.^6 - 4 * rho.^4) .* sin(3 * theta);
                case 17
                    % Z17: Secondary astigmatism 0 deg (higher)
                    Z = sqrt(14) * (15 * rho.^8 - 20 * rho.^6 + 6 * rho.^4) .* cos(2 * theta);
                case 18
                    % Z18: Secondary astigmatism 45 deg (higher)
                    Z = sqrt(14) * (15 * rho.^8 - 20 * rho.^6 + 6 * rho.^4) .* sin(2 * theta);
                case 19
                    % Z19: Tertiary coma X
                    Z = sqrt(14) * (35 * rho.^8 - 60 * rho.^6 + 30 * rho.^4 - 4 * rho.^2) .* rho .* cos(theta);
                case 20
                    % Z20: Tertiary coma Y
                    Z = sqrt(14) * (35 * rho.^8 - 60 * rho.^6 + 30 * rho.^4 - 4 * rho.^2) .* rho .* sin(theta);
                case 21
                    % Z21: Tertiary spherical
                    Z = sqrt(14) * (70 * rho.^8 - 105 * rho.^6 + 42 * rho.^4 - 6 * rho.^2 + 1);
                case 22
                    % Z22: Quadrafoil X
                    Z = sqrt(14) * (21 * rho.^8 - 14 * rho.^6 + 4 * rho.^4) .* cos(4 * theta);
                case 23
                    % Z23: Quadrafoil Y
                    Z = sqrt(14) * (21 * rho.^8 - 14 * rho.^6 + 4 * rho.^4) .* sin(4 * theta);
                case 24
                    % Z24: Secondary trefoil X
                    Z = sqrt(16) * (56 * rho.^10 - 84 * rho.^8 + 36 * rho.^6 - 4 * rho.^4) .* cos(3 * theta);
                case 25
                    % Z25: Secondary trefoil Y
                    Z = sqrt(16) * (56 * rho.^10 - 84 * rho.^8 + 36 * rho.^6 - 4 * rho.^4) .* sin(3 * theta);
                case 26
                    % Z26: Quaternary astigmatism 0 deg
                    Z = sqrt(16) * (28 * rho.^10 - 42 * rho.^8 + 18 * rho.^6 - 2 * rho.^4) .* cos(2 * theta);
                case 27
                    % Z27: Quaternary astigmatism 45 deg
                    Z = sqrt(16) * (28 * rho.^10 - 42 * rho.^8 + 18 * rho.^6 - 2 * rho.^4) .* sin(2 * theta);
                case 28
                    % Z28: Quaternary coma X
                    Z = sqrt(16) * (126 * rho.^12 - 252 * rho.^10 + 168 * rho.^8 - 48 * rho.^6 + 5 * rho.^4) .* rho .* cos(theta);
                case 29
                    % Z29: Quaternary coma Y
                    Z = sqrt(16) * (126 * rho.^12 - 252 * rho.^10 + 168 * rho.^8 - 48 * rho.^6 + 5 * rho.^4) .* rho .* sin(theta);
                case 30
                    % Z30: Quaternary spherical
                    Z = sqrt(16) * (252 * rho.^12 - 504 * rho.^10 + 378 * rho.^8 - 144 * rho.^6 + 24 * rho.^4 - 1);
                case 31
                    % Z31: Pentafoil X
                    Z = sqrt(16) * (126 * rho.^12 - 180 * rho.^10 + 90 * rho.^8 - 20 * rho.^6 + 5 * rho.^4) .* cos(5 * theta);
                case 32
                    % Z32: Pentafoil Y
                    Z = sqrt(16) * (126 * rho.^12 - 180 * rho.^10 + 90 * rho.^8 - 20 * rho.^6 + 5 * rho.^4) .* sin(5 * theta);
                case 33
                    % Z33: Secondary trefoil X (higher)
                    Z = sqrt(18) * (210 * rho.^12 - 378 * rho.^10 + 210 * rho.^8 - 45 * rho.^6 + 4 * rho.^4) .* cos(3 * theta);
                case 34
                    % Z34: Secondary trefoil Y (higher)
                    Z = sqrt(18) * (210 * rho.^12 - 378 * rho.^10 + 210 * rho.^8 - 45 * rho.^6 + 4 * rho.^4) .* sin(3 * theta);
                case 35
                    % Z35: Secondary quadrafoil X
                    Z = sqrt(18) * (495 * rho.^14 - 990 * rho.^12 + 693 * rho.^10 - 220 * rho.^8 + 30 * rho.^6 - 2 * rho.^4) .* cos(4 * theta);
                case 36
                    % Z36: Secondary quadrafoil Y
                    Z = sqrt(18) * (495 * rho.^14 - 990 * rho.^12 + 693 * rho.^10 - 220 * rho.^8 + 30 * rho.^6 - 2 * rho.^4) .* sin(4 * theta);
                otherwise
                    error('ZernikeUtils:invalidIndex', ...
                        'Zernike index n must be 1-36, got %d', n);
            end
        end

        function name = zernikeName(n)
            % zernikeName - Return the name string for Zernike index n (Noll order)

            names = {
                'Piston', ...
                'Tilt X', ...
                'Tilt Y', ...
                'Defocus', ...
                'Astigmatism 0°', ...
                'Astigmatism 45°', ...
                'Coma X', ...
                'Coma Y', ...
                'Spherical', ...
                'Secondary Astigmatism 0°', ...
                'Secondary Astigmatism 45°', ...
                'Secondary Coma X', ...
                'Secondary Coma Y', ...
                'Secondary Spherical', ...
                'Trefoil X', ...
                'Trefoil Y', ...
                'Secondary Astigmatism 0° (high)', ...
                'Secondary Astigmatism 45° (high)', ...
                'Tertiary Coma X', ...
                'Tertiary Coma Y', ...
                'Tertiary Spherical', ...
                'Quadrafoil X', ...
                'Quadrafoil Y', ...
                'Secondary Trefoil X', ...
                'Secondary Trefoil Y', ...
                'Quaternary Astigmatism 0°', ...
                'Quaternary Astigmatism 45°', ...
                'Quaternary Coma X', ...
                'Quaternary Coma Y', ...
                'Quaternary Spherical', ...
                'Pentafoil X', ...
                'Pentafoil Y', ...
                'Secondary Trefoil X (high)', ...
                'Secondary Trefoil Y (high)', ...
                'Secondary Quadrafoil X', ...
                'Secondary Quadrafoil Y'
            };

            if n < 1 || n > 36
                error('ZernikeUtils:invalidIndex', ...
                    'Zernike index n must be 1-36, got %d', n);
            end

            name = names{n};
        end

        function M = zernikeMatrix(rho, theta, nTerms)
            % zernikeMatrix - Build design matrix for Zernike fitting
            %
            % Parameters:
            %   rho    - [Ny x Nx] radial coordinates (normalized 0-1)
            %   theta  - [Ny x Nx] azimuthal coordinates (radians)
            %   nTerms - number of Zernike terms to include (1-36)
            %
            % Returns:
            %   M      - [Ny*Nx x nTerms] design matrix where M(:, n) = Z_n(rho, theta)

            if nTerms < 1 || nTerms > 36
                error('ZernikeUtils:invalidNTerms', ...
                    'nTerms must be 1-36, got %d', nTerms);
            end

            % Reshape to column vectors for matrix construction
            rho_col = rho(:);
            theta_col = theta(:);
            nPixels = length(rho_col);

            % Build design matrix
            M = zeros(nPixels, nTerms);
            for n = 1:nTerms
                M(:, n) = ZernikeUtils.zernike(n, rho_col, theta_col);
            end
        end
    end
end