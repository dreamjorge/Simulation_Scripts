classdef BeamSimulation < handle
%BEAMSIMULATION Base class for beam propagation simulations
%   This class provides common functionality for beam simulation scripts,
%   including parameter management, grid creation, and propagation methods.
%
%   Example:
%       sim = BeamSimulation;
%       sim.setPhysicalParameters('w0', 100e-6, 'lambda', 632.8e-9);
%       sim.createGrid();
%       field = sim.computeInitialField();
    
    properties
        % Physical parameters
        InitialWaist           % Initial beam waist (m)
        Wavelength             % Wavelength (m)
        zCoordinate            % Current z position (m)
        
        % Simulation parameters
        Nx                     % Number of points in x
        Ny                     % Number of points in y  
        Nz                     % Number of points in z
        Dx                     % Window size x (m)
        Dy                     % Window size y (m)
        Dz                     % Propagation window (m)
        
        % Grid arrays (computed)
        x                      % x spatial vector
        y                      % y spatial vector
        z                      % z spatial vector
        kx                     % x wavenumber vector
        ky                     % y wavenumber vector
        
        % Objects
        Grid                   % GridUtils instance
        FFT                    % FFTUtils instance
        
        % Visualization
        Colormap               % Custom colormap
        GenerateVideo         % Video generation flag
    end
    
    methods
        function self = BeamSimulation()
            %BEAMSIMULATION Constructor
            
            % Default physical parameters (HeNe laser)
            self.InitialWaist = 100e-6;  % 100 microns
            self.Wavelength = 632.8e-9; % 632.8 nm
            
            % Default simulation parameters
            self.Nx = 2^10;
            self.Ny = 2^10;
            self.Nz = 2^7;
            
            % Initialize utilities
            self.Grid = GridUtils(self.Nx, self.Nx, 1, 1);
            self.FFT = FFTUtils(true, true);
            
            % Default visualization
            self.Colormap = PhysicalConstants;  % Placeholder
            self.GenerateVideo = false;
        end
        
        function setPhysicalParameters(self, varargin)
            %SETPHYSICALPARAMETERS Set physical parameters
            %   setPhysicalParameters('w0', value, 'lambda', value)
            %
            %   Parameters:
            %     'w0'     - Initial waist (m)
            %     'lambda' - Wavelength (m)
            %     'z'      - z coordinate (m)
            
            for i = 1:2:nargin-1
                param = varargin{i};
                value = varargin{i+1};
                
                switch param
                    case 'w0'
                        self.InitialWaist = value;
                    case 'lambda'
                        self.Wavelength = value;
                    case 'z'
                        self.zCoordinate = value;
                    otherwise
                        warning('BeamSimulation:UnknownParam', ...
                            'Unknown parameter: %s', param);
                end
            end
        end
        
        function createGrid(self, varargin)
            %CREATEGRID Create computational grid
            %   createGrid() uses default parameters
            %   createGrid('Nx', 1024, 'Dz', 0.1) customizes
            
            % Parse optional arguments
            if nargin >= 2
                for i = 1:2:nargin-1
                    param = varargin{i};
                    value = varargin{i+1};
                    switch param
                        case 'Nx', self.Nx = value;
                        case 'Ny', self.Ny = value;
                        case 'Nz', self.Nz = value;
                        case 'Dx', self.Dx = value;
                        case 'Dy', self.Dy = value;
                        case 'Dz', self.Dz = value;
                    end
                end
            end
            
            % Calculate derived parameters
            zr = PhysicalConstants.rayleighDistance(self.InitialWaist, self.Wavelength);
            
            if isempty(self.Dz)
                self.Dz = zr;
            end
            if isempty(self.Dx)
                % Estimate max waist in propagation window
                maxWaist = PhysicalConstants.waistAtZ(self.InitialWaist, self.Dz, ...
                    self.Wavelength, zr);
                self.Dx = 1.2 * 2 * maxWaist;
            end
            self.Dy = self.Dx;
            
            % Update GridUtils
            self.Grid = GridUtils(self.Nx, self.Nx, self.Dx, self.Dx, self.Nz, self.Dz);
            
            % Create spatial vectors
            n = -self.Nx/2:self.Nx/2-1;
            dx = self.Dx / self.Nx;
            self.x = n * dx;
            self.y = self.x;
            self.z = linspace(0, self.Dz, self.Nz);
            
            % Create frequency vectors
            du = 1 / self.Dx;
            u = n * du;
            self.kx = 2*pi*u;
            self.ky = self.kx;
        end
        
        function info = getInfo(self)
            %GETINFO Return simulation info as struct
            
            zr = PhysicalConstants.rayleighDistance(self.InitialWaist, self.Wavelength);
            
            info = struct(...
                'InitialWaist', self.InitialWaist, ...
                'Wavelength', self.Wavelength, ...
                'RayleighDistance', zr, ...
                'Nx', self.Nx, ...
                'Ny', self.Ny, ...
                'Nz', self.Nz, ...
                'Dx', self.Dx, ...
                'Dz', self.Dz);
        end
        
        function validate(self)
            %VALIDATE Check simulation parameters
            
            if isempty(self.InitialWaist) || self.InitialWaist <= 0
                error('BeamSimulation:InvalidWaist', ...
                    'InitialWaist must be positive');
            end
            
            if isempty(self.Wavelength) || self.Wavelength <= 0
                error('BeamSimulation:InvalidWavelength', ...
                    'Wavelength must be positive');
            end
            
            if self.Nx < 64 || self.Ny < 64
                warning('BeamSimulation:LowResolution', ...
                    'N < 64 may cause sampling artifacts');
            end
        end
    end
    
    methods (Abstract)
        % Subclasses must implement these
        field = computeInitialField(self, varargin)
        field = propagate(self, z)
    end
end
