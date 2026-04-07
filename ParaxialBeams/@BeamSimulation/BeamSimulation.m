classdef BeamSimulation < handle
%BEAMSIMULATION Base class for beam propagation simulations
    
    properties
        InitialWaist, Wavelength, zCoordinate
        Nx=2^10, Ny=2^10, Nz=2^7
        Dx, Dy, Dz
        x, y, z, kx, ky
        Grid, FFT
    end
    
    methods
        function self = BeamSimulation()
            self.InitialWaist = 100e-6;
            self.Wavelength = 632.8e-9;
            self.Grid = GridUtils(self.Nx, self.Nx, 1, 1);
            self.FFT = FFTUtils(true, true);
        end
        
        function setPhysicalParameters(self, varargin)
            for i = 1:2:nargin-1
                param = varargin{i}; value = varargin{i+1};
                switch param
                    case 'w0', self.InitialWaist = value;
                    case 'lambda', self.Wavelength = value;
                    case 'z', self.zCoordinate = value;
                end
            end
        end
        
        function createGrid(self, varargin)
            zr = PhysicalConstants.rayleighDistance(self.InitialWaist, self.Wavelength);
            if isempty(self.Dz), self.Dz = zr; end
            maxWaist = PhysicalConstants.waistAtZ(self.InitialWaist, self.Dz, self.Wavelength, zr);
            if isempty(self.Dx), self.Dx = 1.2 * 2 * maxWaist; end
            self.Dy = self.Dx;
            
            self.Grid = GridUtils(self.Nx, self.Nx, self.Dx, self.Dx, self.Nz, self.Dz);
            
            n = -self.Nx/2:self.Nx/2-1;
            dx = self.Dx / self.Nx;
            self.x = n * dx; self.y = self.x;
            self.z = linspace(0, self.Dz, self.Nz);
            
            du = 1 / self.Dx;
            u = n * du;
            self.kx = 2*pi*u; self.ky = self.kx;
        end
    end
end
