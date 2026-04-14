% Compatible with GNU Octave and MATLAB
% Tests for OpticalRay

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== OpticalRay Tests ===\n\n');
passed = 0;
failed = 0;

% testDefaultConstructor
try
    ray = OpticalRay();
    if isscalar(ray)
        fprintf('  PASS: default constructor\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: default constructor not scalar\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testDefaultProperties
try
    ray = OpticalRay();
    props = {'xCoordinate', 'yCoordinate', 'zCoordinate', ...
            'zxSlope', 'zySlope', 'xySlope', 'hankelType'};
    all_have = true;
    for i = 1:numel(props)
        if ~isprop(ray, props{i})
            all_have = false;
            break;
        end
    end
    if all_have
        fprintf('  PASS: all properties exist\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: missing properties\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testDefaultCoordinateValues
try
    ray = OpticalRay();
    if ray.xCoordinate == 0 && ray.yCoordinate == 0 && ray.zCoordinate == 0
        fprintf('  PASS: default coordinate values\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: default coordinate values incorrect\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testDefaultSlopes
try
    ray = OpticalRay();
    if ray.zxSlope == Inf && ray.zySlope == Inf && ray.xySlope == Inf
        fprintf('  PASS: default slopes are Inf\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: default slopes should be Inf\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testDefaultHankelType
try
    ray = OpticalRay();
    if ray.hankelType == 1
        fprintf('  PASS: default hankelType is 1\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: default hankelType should be 1\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testSetCoordinates
try
    ray = OpticalRay();
    ray.xCoordinate = 0.5;
    ray.yCoordinate = 0.3;
    ray.zCoordinate = 0.1;
    if ray.xCoordinate == 0.5 && ray.yCoordinate == 0.3 && ray.zCoordinate == 0.1
        fprintf('  PASS: set coordinates\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: set coordinates\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testSetSlopes
try
    ray = OpticalRay();
    ray.zxSlope = 10;
    ray.zySlope = 5;
    ray.xySlope = 2;
    if ray.zxSlope == 10 && ray.zySlope == 5 && ray.xySlope == 2
        fprintf('  PASS: set slopes\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: set slopes\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testSetHankelType
try
    ray = OpticalRay();
    ray.hankelType = 2;
    if ray.hankelType == 2
        fprintf('  PASS: set hankelType\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: set hankelType\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

fprintf('\n=== OpticalRay: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end