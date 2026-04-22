function config = bench_config(mode)
% bench_config - Benchmark matrix/configuration for propagator performance.
% Compatible with GNU Octave and MATLAB.

    if nargin < 1 || isempty(mode)
        mode = 'quick';
    end

    mode = lower(mode);
    if ~ismember(mode, {'quick', 'full'})
        error('bench_config:invalidMode', 'Mode must be ''quick'' or ''full''.');
    end

    config = struct();
    config.mode = mode;
    config.seed = 12345;
    config.w0 = 100e-6;
    config.lambda = 632.8e-9;
    config.outputPrefix = 'baseline';
    config.outputDir = fullfile('docs', 'performance');

    config.grid.small = struct('tier', 'small', 'Nx', 256, 'Ny', 256, 'Dx', 1e-3, 'Dy', 1e-3);
    config.grid.medium = struct('tier', 'medium', 'Nx', 512, 'Ny', 512, 'Dx', 1e-3, 'Dy', 1e-3);
    config.grid.large = struct('tier', 'large', 'Nx', 1024, 'Ny', 1024, 'Dx', 1e-3, 'Dy', 1e-3);

    if strcmp(mode, 'quick')
        config.warmup = 1;
        config.repeats = 3;
        config.scenarios = {
            struct('name', 'gaussian_fft_small',      'propagator', 'fft',      'beamType', 'gaussian', 'beamParams', struct(),              'gridTier', 'small',  'zFinal', 0.05, 'rayDz', []), ...
            struct('name', 'gaussian_analytic_small', 'propagator', 'analytic', 'beamType', 'gaussian', 'beamParams', struct(),              'gridTier', 'small',  'zFinal', 0.05, 'rayDz', []), ...
            struct('name', 'gaussian_raytrace_small', 'propagator', 'raytrace', 'beamType', 'gaussian', 'beamParams', struct(),              'gridTier', 'small',  'zFinal', 0.01, 'rayDz', 1e-4)
        };
    else
        config.warmup = 2;
        config.repeats = 5;
        config.scenarios = {
            struct('name', 'gaussian_fft_small',        'propagator', 'fft',      'beamType', 'gaussian', 'beamParams', struct(),                   'gridTier', 'small',  'zFinal', 0.05, 'rayDz', []), ...
            struct('name', 'gaussian_fft_medium',       'propagator', 'fft',      'beamType', 'gaussian', 'beamParams', struct(),                   'gridTier', 'medium', 'zFinal', 0.05, 'rayDz', []), ...
            struct('name', 'gaussian_fft_large',        'propagator', 'fft',      'beamType', 'gaussian', 'beamParams', struct(),                   'gridTier', 'large',  'zFinal', 0.05, 'rayDz', []), ...
            struct('name', 'gaussian_analytic_small',   'propagator', 'analytic', 'beamType', 'gaussian', 'beamParams', struct(),                   'gridTier', 'small',  'zFinal', 0.05, 'rayDz', []), ...
            struct('name', 'gaussian_analytic_medium',  'propagator', 'analytic', 'beamType', 'gaussian', 'beamParams', struct(),                   'gridTier', 'medium', 'zFinal', 0.05, 'rayDz', []), ...
            struct('name', 'gaussian_analytic_large',   'propagator', 'analytic', 'beamType', 'gaussian', 'beamParams', struct(),                   'gridTier', 'large',  'zFinal', 0.05, 'rayDz', []), ...
            struct('name', 'gaussian_raytrace_small',   'propagator', 'raytrace', 'beamType', 'gaussian', 'beamParams', struct(),                   'gridTier', 'small',  'zFinal', 0.01, 'rayDz', 1e-4), ...
            struct('name', 'gaussian_raytrace_medium',  'propagator', 'raytrace', 'beamType', 'gaussian', 'beamParams', struct(),                   'gridTier', 'medium', 'zFinal', 0.01, 'rayDz', 1e-4), ...
            struct('name', 'hermite_fft_medium',        'propagator', 'fft',      'beamType', 'hermite',  'beamParams', struct('n', 1, 'm', 1),    'gridTier', 'medium', 'zFinal', 0.05, 'rayDz', []), ...
            struct('name', 'hermite_analytic_medium',   'propagator', 'analytic', 'beamType', 'hermite',  'beamParams', struct('n', 1, 'm', 1),    'gridTier', 'medium', 'zFinal', 0.05, 'rayDz', []), ...
            struct('name', 'laguerre_fft_medium',       'propagator', 'fft',      'beamType', 'laguerre', 'beamParams', struct('l', 1, 'p', 0),    'gridTier', 'medium', 'zFinal', 0.05, 'rayDz', []), ...
            struct('name', 'laguerre_analytic_medium',  'propagator', 'analytic', 'beamType', 'laguerre', 'beamParams', struct('l', 1, 'p', 0),    'gridTier', 'medium', 'zFinal', 0.05, 'rayDz', [])
        };
    end
end
