classdef VisualizationUtils
    % VisualizationUtils Modernized visualization tools for ray tracing.
    
    methods (Static)
        function plotRays3D(bundle, varargin)
            % Plot rays in 3D
            figure;
            hold on;
            for i = 1:bundle.Ny
                for j = 1:bundle.Nx
                    plot3(squeeze(bundle.z(i,j,:)), squeeze(bundle.x(i,j,:)), squeeze(bundle.y(i,j,:)), varargin{:});
                end
            end
            grid on;
            view(3);
            xlabel('z (m)'); ylabel('x (m)'); zlabel('y (m)');
            title('Ray Tracing Propagation (3D)');
        end
        
        function plotRays2D(bundle, plane, varargin)
            % Plot rays in 2D (xz or yz plane)
            if nargin < 2, plane = 'xz'; end
            figure;
            hold on;
            if strcmpi(plane, 'xz')
                for i = 1:bundle.Ny
                    for j = 1:bundle.Nx
                        plot(squeeze(bundle.z(i,j,:)), squeeze(bundle.x(i,j,:)), varargin{:});
                    end
                end
                ylabel('x (m)');
            else
                for i = 1:bundle.Ny
                    for j = 1:bundle.Nx
                        plot(squeeze(bundle.z(i,j,:)), squeeze(bundle.y(i,j,:)), varargin{:});
                    end
                end
                ylabel('y (m)');
            end
            grid on;
            xlabel('z (m)');
            title(sprintf('Ray Tracing Propagation (%s)', upper(plane)));
        end
    end
end
