function uninstall()
    % uninstall - Uninstall Simulation_Scripts package
    %
    % This script is called automatically by:
    %   Octave: pkg uninstall simulation_scripts
    %   MATLAB: matlab.addons.uninstall('Simulation_Scripts')
    %
    % Or manually:
    %   octave --eval "uninstall"
    %   matlab >> run uninstall.m

    scriptPath = fileparts(mfilename('fullpath'));

    % Remove main repo root
    rmpath(scriptPath);

    % Remove all +paraxial/ namespace subpackages
    paraxialRoot = fullfile(scriptPath, '+paraxial');
    subfolders = {...
        '+beams', ...
        '+parameters', ...
        '+computation', ...
        '+propagation', ...
        '+propagation/+field', ...
        '+propagation/+rays', ...
        '+visualization'};
    for i = 1:numel(subfolders)
        subPath = fullfile(paraxialRoot, subfolders{i});
        if exist(subPath, 'dir')
            rmpath(subPath);
        end
    end

    % Remove the paraxial root last
    if exist(paraxialRoot, 'dir')
        rmpath(paraxialRoot);
    end

    fprintf('[Simulation_Scripts] Uninstalled successfully.\n');
end