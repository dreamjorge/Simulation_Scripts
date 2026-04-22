function summary = run_benchmarks(mode, outputDir)
% run_benchmarks - Deterministic benchmark runner for core propagators.
% Compatible with GNU Octave and MATLAB.
%
% Usage:
%   summary = run_benchmarks();
%   summary = run_benchmarks('quick');
%   summary = run_benchmarks('full');
%   summary = run_benchmarks('full', fullfile('docs', 'performance'));

    if nargin < 1 || isempty(mode)
        mode = 'quick';
    end

    cfg = bench_config(mode);
    if nargin >= 2 && ~isempty(outputDir)
        cfg.outputDir = outputDir;
    end

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(repoRoot, 'src', 'beams'));
    addpath(fullfile(repoRoot, 'src', 'parameters'));
    addpath(fullfile(repoRoot, 'src', 'computation'));
    addpath(fullfile(repoRoot, 'src', 'propagation', 'field'));
    addpath(fullfile(repoRoot, 'src', 'propagation', 'rays'));
    addpath(fullfile(repoRoot, 'src', 'visualization'));
    addpath(fullfile(repoRoot, 'ParaxialBeams'));
    addpath(fullfile(repoRoot, 'ParaxialBeams', 'Addons'));
    addpath(fileparts(mfilename('fullpath')));

    utils = bench_utils();
    utils.setSeed(cfg.seed);

    absOutputDir = fullfile(repoRoot, cfg.outputDir);
    utils.ensureDir(absOutputDir);

    fprintf('=== Propagator Benchmark Runner (%s mode) ===\n', upper(cfg.mode));
    fprintf('Scenarios: %d | Warm-up: %d | Repeats: %d\n\n', ...
        numel(cfg.scenarios), cfg.warmup, cfg.repeats);

    rows = cell(numel(cfg.scenarios), 1);
    for i = 1:numel(cfg.scenarios)
        sc = cfg.scenarios{i};
        gridSpec = cfg.grid.(sc.gridTier);
        gridObj = GridUtils(gridSpec.Nx, gridSpec.Ny, gridSpec.Dx, gridSpec.Dy);
        beamArgs = name_value_args(sc.beamParams);
        beamObj = BeamFactory.create(sc.beamType, cfg.w0, cfg.lambda, beamArgs{:});

        runFcn = @() run_single_scenario(sc, gridObj, beamObj, cfg.lambda);
        [medianSec, samples] = utils.timeScenario(runFcn, cfg.warmup, cfg.repeats);

        row = struct();
        row.timestamp = utils.timestamp();
        row.mode = cfg.mode;
        row.scenario = sc.name;
        row.propagator = sc.propagator;
        row.beamType = sc.beamType;
        row.gridTier = sc.gridTier;
        row.nx = gridSpec.Nx;
        row.ny = gridSpec.Ny;
        row.zFinal = sc.zFinal;
        row.warmup = cfg.warmup;
        row.repeats = cfg.repeats;
        row.medianRuntimeSec = medianSec;
        row.minRuntimeSec = min(samples);
        row.maxRuntimeSec = max(samples);
        rows{i} = row;

        fprintf('[%2d/%2d] %-30s median = %.6f s\n', ...
            i, numel(cfg.scenarios), sc.name, medianSec);
    end

    runStamp = datestr(now, 'yyyymmdd_HHMMSS');
    csvName = sprintf('%s_%s_%s.csv', cfg.outputPrefix, cfg.mode, runStamp);
    csvPath = fullfile(absOutputDir, csvName);
    utils.writeCsv(csvPath, rows);

    if strcmp(cfg.mode, 'full')
        baselineName = sprintf('baseline_%s.csv', datestr(now, 'yyyy-mm-dd'));
        baselinePath = fullfile(absOutputDir, baselineName);
        utils.writeCsv(baselinePath, rows);
    else
        baselinePath = '';
    end

    fprintf('\nCSV written: %s\n', csvPath);
    if ~isempty(baselinePath)
        fprintf('Baseline written: %s\n', baselinePath);
    end

    summary = struct();
    summary.mode = cfg.mode;
    summary.scenarioCount = numel(rows);
    summary.csvPath = csvPath;
    summary.baselinePath = baselinePath;
    summary.rows = rows;
end

function args = name_value_args(s)
    if isempty(s)
        args = {};
        return;
    end

    fn = fieldnames(s);
    args = cell(1, numel(fn) * 2);
    idx = 1;
    for i = 1:numel(fn)
        args{idx} = fn{i};
        args{idx + 1} = s.(fn{i});
        idx = idx + 2;
    end
end

function run_single_scenario(sc, gridObj, beamObj, lambda)
    switch lower(sc.propagator)
        case 'analytic'
            prop = AnalyticPropagator(gridObj);
            field = prop.propagate(beamObj, sc.zFinal); %#ok<NASGU>

        case 'fft'
            prop = FFTPropagator(gridObj, lambda);
            field = prop.propagate(beamObj, sc.zFinal); %#ok<NASGU>

        case 'raytrace'
            if isempty(sc.rayDz)
                prop = RayTracePropagator(gridObj, 'RK4');
            else
                prop = RayTracePropagator(gridObj, 'RK4', sc.rayDz);
            end
            bundle = prop.propagate(beamObj, sc.zFinal); %#ok<NASGU>

        otherwise
            error('run_benchmarks:unknownPropagator', ...
                'Unknown propagator type: %s', sc.propagator);
    end
end
