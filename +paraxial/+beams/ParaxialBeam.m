classdef (Abstract) ParaxialBeam < handle
    % ParaxialBeam - Abstract base class for all paraxial beam models.
    %
    % Every beam type MUST inherit from this class and implement the three
    % interface methods below. This enforces a unified API that lets
    % propagators (FFT, ray tracing, analytic) work with any beam without
    % knowing its internal coordinate system or formula.
    %
    % Octave/MATLAB portability note:
    %   MATLAB supports 'methods (Abstract)' with bodyless signatures.
    %   Octave requires method bodies even for abstract contracts, so the
    %   three interface methods below throw an explicit error. The class-level
    %   (Abstract) attribute still prevents direct instantiation.
    %
    % Interface (methods every subclass must override):
    %
    %   field = opticalField(obj, X, Y, z)
    %     Returns the complex optical field on a Cartesian grid (X, Y) at
    %     propagation distance z. X and Y are 2-D matrices in metres.
    %     Subclasses that work in polar coordinates convert internally.
    %
    %   params = getParameters(obj, z)
    %     Returns a GaussianParameters object evaluated at axial position z.
    %     Useful for propagators that need beam radius, Gouy phase, etc.
    %
    %   name = beamName(obj)
    %     Returns a short string identifying the beam type, e.g. 'gaussian',
    %     'hermite_3_2', 'laguerre_2_1'. Used by BeamFactory and logging.
    %
    % Shared state (set by this constructor, available to all subclasses):
    %   Lambda  - wavelength (m)
    %   k       - wave number 2*pi/lambda (rad/m)

    properties
        Lambda  % Wavelength (m)
        k       % Wave number (rad/m)
    end

    methods
        function obj = ParaxialBeam(lambda)
            if nargin > 0
                obj.Lambda = lambda;
                obj.k      = 2 * pi / lambda;
            end
        end
    end

    % These stubs enforce the interface contract. Subclasses MUST override
    % all three. The error message names the missing method and the class.
    methods
        function field = opticalField(obj, X, Y, z)
            % opticalField - Complex field on a Cartesian (X,Y) grid at depth z.
            % Subclasses must override this method.
            error('ParaxialBeam:notImplemented', ...
                '%s must implement opticalField(obj, X, Y, z).', class(obj));
        end

        function params = getParameters(obj, z)
            % getParameters - GaussianParameters evaluated at axial position z.
            % Subclasses must override this method.
            error('ParaxialBeam:notImplemented', ...
                '%s must implement getParameters(obj, z).', class(obj));
        end

        function name = beamName(obj)
            % beamName - Short string identifier for this beam type.
            % Subclasses must override this method.
            error('ParaxialBeam:notImplemented', ...
                '%s must implement beamName(obj).', class(obj));
        end
    end
end
