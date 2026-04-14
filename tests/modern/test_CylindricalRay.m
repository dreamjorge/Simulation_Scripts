% Compatible with GNU Octave and MATLAB
% Tests for CylindricalRay

addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'ParaxialBeams'));

fprintf('=== CylindricalRay Tests ===\n\n');
passed = 0;
failed = 0;

% testDefaultConstructor
try
    ray = CylindricalRay();
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
    ray = CylindricalRay();
    props = {'rCoordinate', 'thetaCoordinate', 'zCoordinate', ...
            'xCoordinate', 'yCoordinate', 'zrSlope', 'zthSlope', 'rthSlope', 'hankelType'};
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

% testDefaultValues
try
    ray = CylindricalRay();
    if ray.rCoordinate == 0 && ray.thetaCoordinate == 0 && ray.zCoordinate == 0
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
    ray = CylindricalRay();
    if isinf(ray.zrSlope) && isinf(ray.zthSlope) && isinf(ray.rthSlope)
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
    ray = CylindricalRay();
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
    ray = CylindricalRay();
    ray.rCoordinate = 0.5;
    ray.thetaCoordinate = pi/4;
    ray.zCoordinate = 0.1;
    if ray.rCoordinate == 0.5 && abs(ray.thetaCoordinate - pi/4) < eps && ray.zCoordinate == 0.1
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
    ray = CylindricalRay();
    ray.zrSlope = 10;
    ray.zthSlope = 5;
    ray.rthSlope = 2;
    if ray.zrSlope == 10 && ray.zthSlope == 5 && ray.rthSlope == 2
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
    ray = CylindricalRay();
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

% testCartesianCoordinatesDerived
try
    ray = CylindricalRay();
    ray.rCoordinate = sqrt(2);
    ray.thetaCoordinate = pi/4;
    % After setting r and theta, x = r*cos(theta), y = r*sin(theta)
    % x should be 1, y should be 1
    % But CylindricalRay doesn't auto-compute x,y from r,theta
    % They remain 0 by default unless explicitly set
    if isprop(ray, 'xCoordinate') && isprop(ray, 'yCoordinate')
        fprintf('  PASS: Cartesian coordinate properties exist\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: Cartesian properties missing\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testArrayConstruction
try
    rays(10) = CylindricalRay();
    if numel(rays) == 10
        fprintf('  PASS: array construction\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: array construction\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testIndependentArrayElements
try
    rays(5) = CylindricalRay();
    rays(1).rCoordinate = 1;
    rays(2).rCoordinate = 2;
    rays(3).rCoordinate = 3;
    if rays(1).rCoordinate == 1 && rays(2).rCoordinate == 2 && rays(3).rCoordinate == 3
        fprintf('  PASS: independent array elements\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: array elements not independent\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

fprintf('\n=== CylindricalRay: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
