% utils/VisualizationUtils.m
function plotBeamIntensity(beam, z, grid_size)
    [x, y] = meshgrid(linspace(grid_size(1), grid_size(2), 100), ...
                      linspace(grid_size(1), grid_size(2), 100));
    
    % Calculate the intensity profile
    intensity = beam.intensity_profile(z, x, y);
    
    % Plot the intensity profile
    figure;
    surf(x, y, intensity);
    shading interp;
    colormap(hot);
    title(sprintf('Intensity Profile at z = %.2f m', z));
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Intensity');
end
