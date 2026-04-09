classdef test_RayTracing < matlab.unittest.TestCase
    % Unit tests for Ray Tracing integration
    
    methods (Test)
        function testRayBundleInitialization(testCase)
            % Test grid creation
            bundle = RayBundle.createGrid(10, 12, 1e-3, 1.2e-3);
            testCase.verifySize(bundle.x, [12, 10]);
            testCase.verifyEqual(max(bundle.x(:)), 0.5e-3, 'RelTol', 1e-5);
            
            % Test concentric creation
            bundle = RayBundle.createConcentric(5, 8, 1e-3);
            testCase.verifyEqual(bundle.Nx, 8);
            testCase.verifyEqual(bundle.Ny, 5);
        end
        
        function testRayTracerSlopes(testCase)
            % Test slope calculation for a Gaussian beam
            lambda = 632.8e-9;
            w0 = 100e-6;
            beam = GaussianBeam(w0, lambda);
            
            % At z=0, the phase of a Gaussian beam is flat (m=0)
            [sx, sy] = RayTracer.calculateSlopes(beam, 0, 0, 0);
            testCase.verifyEqual(sx, 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(sy, 0, 'AbsTol', 1e-10);
            
            % Off-center at z=0 should also be 0 in paraxial approx if tilted=false
            [sx, sy] = RayTracer.calculateSlopes(beam, 50e-6, 0, 0);
            testCase.verifyEqual(sx, 0, 'AbsTol', 1e-6);
        end
        
        function testPropagationFreeSpace(testCase)
            % Test Euler and RK4 for a Gaussian beam
            lambda = 632.8e-9;
            w0 = 100e-6;
            beam = GaussianBeam(w0, lambda);
            zr = pi * w0^2 / lambda;
            
            bundle = RayBundle.createGrid(5, 5, 200e-6, 200e-6);
            
            % Propagate to Rayleigh distance
            dz = zr / 10;
            RayTracer.propagate(bundle, beam, zr, dz, 'RK4');
            
            % Rays should diverge
            testCase.verifyGreaterThan(bundle.Nz, 1);
            last_x = bundle.x(1, end, end);
            first_x = bundle.x(1, end, 1);
            testCase.verifyGreaterThan(abs(last_x), abs(first_x));
        end
    end
end
