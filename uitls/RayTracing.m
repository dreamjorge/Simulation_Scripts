% utils/RayTracing.m
function ray_tracing(beam, z_values, grid_size)
    % Ray tracing utility for wavefront visualization
    
    [x, y] = meshgrid(linspace(grid_size(1), grid_size(2), 100), ...
                      linspace(grid_size(1), grid_size(2), 100));
    
    figure;
    hold on;
    
    for z = z_values
        % Calculate the radius of curvature and beam waist at each z
        R = beam.radius_of_curvature(z);
        w = beam.beam_waist(z);
        
        % Visualize wavefront by plotting phase profile
        phase_front = sqrt(x.^2 + y.^2) / R;
        surf(x, y, phase_front);
        shading interp;
        colormap(jet);
        title(sprintf('Wavefront at z = %.2f m', z));
        xlabel('x (m)');
        ylabel('y (m)');
        zlabel('Phase');
        pause(0.5);
    end
    
    hold off;
end
