function utils = bench_utils()
% bench_utils - Helper functions used by performance benchmark runner.
% Compatible with GNU Octave and MATLAB.

    utils.setSeed = @setSeed;
    utils.ensureDir = @ensureDir;
    utils.timestamp = @timestamp;
    utils.timeScenario = @timeScenario;
    utils.writeCsv = @writeCsv;
end

function setSeed(seed)
    if nargin < 1 || isempty(seed)
        seed = 12345;
    end

    if exist('rng', 'file') == 2
        rng(seed, 'twister');
    else
        rand('state', seed); %#ok<RAND>
        randn('state', seed); %#ok<RAND>
    end
end

function ensureDir(pathStr)
    if ~exist(pathStr, 'dir')
        mkdir(pathStr);
    end
end

function ts = timestamp()
    ts = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
end

function [medianSec, samplesSec] = timeScenario(runFcn, warmup, repeats)
    if nargin < 2 || isempty(warmup), warmup = 1; end
    if nargin < 3 || isempty(repeats), repeats = 3; end

    for i = 1:warmup
        runFcn();
    end

    samplesSec = zeros(repeats, 1);
    for i = 1:repeats
        t0 = tic;
        runFcn();
        samplesSec(i) = toc(t0);
    end

    medianSec = median(samplesSec);
end

function writeCsv(filePath, rows)
    fid = fopen(filePath, 'w');
    if fid == -1
        error('bench_utils:csvOpenFailed', 'Could not open output file: %s', filePath);
    end

    try
        fprintf(fid, ['timestamp,mode,scenario,propagator,beam_type,grid_tier,' ...
                      'nx,ny,z_final,warmup,repeats,median_runtime_sec,min_runtime_sec,max_runtime_sec\n']);

        for i = 1:numel(rows)
            r = rows{i};
            fprintf(fid, '%s,%s,%s,%s,%s,%s,%d,%d,%.12g,%d,%d,%.12g,%.12g,%.12g\n', ...
                r.timestamp, r.mode, r.scenario, r.propagator, r.beamType, r.gridTier, ...
                r.nx, r.ny, r.zFinal, r.warmup, r.repeats, r.medianRuntimeSec, ...
                r.minRuntimeSec, r.maxRuntimeSec);
        end
    catch ME
        fclose(fid);
        rethrow(ME);
    end

    fclose(fid);
end
