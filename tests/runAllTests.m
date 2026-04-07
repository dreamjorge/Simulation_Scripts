% Test suite for Simulation_Scripts
% Run with: matlab -batch "runtests('tests/')"

%% PhysicalConstants tests
classdef TestPhysicalConstants < matlab.unittest.TestCase
    methods (Test)
        function testWaveNumber(testCase)
            % k = 2*pi/lambda
            lambda = 632.8e-9;
            k = PhysicalConstants.waveNumber(lambda);
            expected = 2*pi / lambda;
            testCase.verifyEqual(k, expected, 'AbsTol', 1e-5);
        end
        
        function testRayleighDistance(testCase)
            % zr = pi*w0^2/lambda
            w0 = 100e-6;
            lambda = 632.8e-9;
            zr = PhysicalConstants.rayleighDistance(w0, lambda);
            expected = pi * w0^2 / lambda;
            testCase.verifyEqual(zr, expected, 'RelTol', 1e-5);
        end
        
        function testWaistAtZ(testCase)
            w0 = 100e-6;
            z = 0.05;
            lambda = 632.8e-9;
            zr = PhysicalConstants.rayleighDistance(w0, lambda);
            w = PhysicalConstants.waistAtZ(w0, z, lambda, zr);
            expected = w0 * sqrt(1 + (z/zr)^2);
            testCase.verifyEqual(w, expected, 'RelTol', 1e-5);
        end
        
        function testWaistAtZEqualsW0AtOrigin(testCase)
            w0 = 100e-6;
            lambda = 632.8e-9;
            w = PhysicalConstants.waistAtZ(w0, 0, lambda);
            testCase.verifyEqual(w, w0, 'RelTol', 1e-5);
        end
    end
end

%% GridUtils tests
classdef TestGridUtils < matlab.unittest.TestCase
    methods (Test)
        function testCreate2DGridSize(testCase)
            Nx = 256; Ny = 128;
            Dx = 1e-3; Dy = 0.5e-3;
            grid = GridUtils(Nx, Ny, Dx, Dy);
            [X, Y] = grid.create2DGrid();
            testCase.verifySize(X, [Ny, Nx]);
            testCase.verifySize(Y, [Ny, Nx]);
        end
        
        function testCreate2DGridValues(testCase)
            Nx = 256; Dx = 1e-3;
            grid = GridUtils(Nx, Nx, Dx, Dx);
            [X, Y] = grid.create2DGrid();
            dx = Dx / Nx;
            % Check center is near zero
            testCase.verifyEqual(X(Nx/2+1, Nx/2+1), 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(Y(Nx/2+1, Nx/2+1), 0, 'AbsTol', 1e-10);
        end
    end
end

%% FFTUtils tests
classdef TestFFTUtils < matlab.unittest.TestCase
    methods (Test)
        function testFFTRoundtrip(testCase)
            fftOps = FFTUtils(true, true);
            % Create test signal
            [X, Y] = meshgrid(linspace(-1,1,64), linspace(-1,1,64));
            R = sqrt(X.^2 + Y.^2);
            g = exp(-R.^2); % Gaussian
            
            g_rec = fftOps.ifft2(fftOps.fft2(g));
            testCase.verifyEqual(g, g_rec, 'AbsTol', 1e-10);
        end
        
        function testFFTNormalized(testCase)
            fftOps = FFTUtils(true, true);
            g = rand(64, 64);
            G = fftOps.fft2(g);
            % Parseval's theorem: sum(|g|^2) = sum(|G|^2)/N
            sum_g = sum(abs(g).^2, 'all');
            sum_G = sum(abs(G).^2, 'all') / numel(g);
            testCase.verifyEqual(sum_g, sum_G, 'RelTol', 1e-10);
        end
        
        function testTransferFunctionAtZero(testCase)
            fftOps = FFTUtils();
            [Kx, Ky] = meshgrid(linspace(-1e6,1e6,32));
            H = fftOps.transferFunction(Kx, Ky, 0, 632.8e-9);
            testCase.verifyEqual(H, ones(32), 'AbsTol', 1e-10);
        end
    end
end

%% GaussianParameters tests
classdef TestGaussianParameters < matlab.unittest.TestCase
    methods (Test)
        function testRayleighDistance(testCase)
            z = 0;
            w0 = 100e-6;
            lambda = 632.8e-9;
            params = GaussianParameters(z, w0, lambda);
            expected_zr = pi * w0^2 / lambda;
            testCase.verifyEqual(params.RayleighDistance, expected_zr, 'RelTol', 1e-5);
        end
        
        function testWaistAtOriginEqualsInitial(testCase)
            z = 0;
            w0 = 100e-6;
            lambda = 632.8e-9;
            params = GaussianParameters(z, w0, lambda);
            testCase.verifyEqual(params.Waist, w0, 'RelTol', 1e-5);
        end
        
        function testWaveNumber(testCase)
            z = 0; w0 = 100e-6; lambda = 632.8e-9;
            params = GaussianParameters(z, w0, lambda);
            expected_k = 2*pi / lambda;
            testCase.verifyEqual(params.k, expected_k, 'RelTol', 1e-5);
        end
        
        function testVectorZCoordinate(testCase)
            z = linspace(0, 0.1, 10);
            w0 = 100e-6;
            lambda = 632.8e-9;
            % Should not error
            params = GaussianParameters(z, w0, lambda);
            testCase.verifyEqual(numel(params.zCoordinate), 10);
        end
    end
end
