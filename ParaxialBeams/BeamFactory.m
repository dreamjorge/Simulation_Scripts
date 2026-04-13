classdef BeamFactory
    % BeamFactory - Factory for creating any ParaxialBeam by name
    % Compatible with GNU Octave and MATLAB
    %
    % Usage:
    %   beam = BeamFactory.create(type, w0, lambda)
    %   beam = BeamFactory.create(type, w0, lambda, Name, Value, ...)
    %
    % Supported types and their optional parameters:
    %
    %   'gaussian'          — GaussianBeam(w0, lambda)
    %   'hermite'           — HermiteBeam(w0, lambda, n, m)
    %                         Name-Value: 'n' (default 0), 'm' (default 0)
    %   'laguerre'          — LaguerreBeam(w0, lambda, l, p)
    %                         Name-Value: 'l' (default 0), 'p' (default 0)
    %   'elegant_hermite'   — ElegantHermiteBeam(w0, lambda, n, m)
    %                         Name-Value: 'n' (default 0), 'm' (default 0)
    %   'elegant_laguerre'  — ElegantLaguerreBeam(w0, lambda, l, p)
    %                         Name-Value: 'l' (default 0), 'p' (default 0)
    %   'hankel'            — HankelLaguerre(w0, lambda, l, p, hankelType)
    %                         Name-Value: 'l' (default 0), 'p' (default 0),
    %                                     'type' (default 1)
    %
    % Examples:
    %   beam = BeamFactory.create('gaussian', 100e-6, 632.8e-9);
    %   beam = BeamFactory.create('hermite', 100e-6, 632.8e-9, 'n', 2, 'm', 1);
    %   beam = BeamFactory.create('laguerre', 100e-6, 632.8e-9, 'l', 1, 'p', 0);
    %   beam = BeamFactory.create('hankel', 100e-6, 632.8e-9, 'l', 2, 'type', 1);
    %
    % Extensibility:
    %   To add a new beam type, add a case to the switch in BeamFactory.create().
    %   No changes to propagators, tests, or calling code required.

    methods (Static)
        function beam = create(type, w0, lambda, varargin)
            % create - Instantiate a ParaxialBeam by name.
            %
            % Parameters:
            %   type   (char):   beam type identifier (see class header)
            %   w0     (scalar): beam waist at z = 0 (m)
            %   lambda (scalar): wavelength (m)
            %   varargin:        Name-Value pairs for mode indices
            %
            % Returns:
            %   beam (ParaxialBeam subclass)

            % Parse optional Name-Value pairs
            n    = BeamFactory.getOpt(varargin, 'n',    0);
            m    = BeamFactory.getOpt(varargin, 'm',    0);
            l    = BeamFactory.getOpt(varargin, 'l',    0);
            p    = BeamFactory.getOpt(varargin, 'p',    0);
            htype = BeamFactory.getOpt(varargin, 'type', 1);

            switch lower(type)
                case 'gaussian'
                    beam = GaussianBeam(w0, lambda);

                case 'hermite'
                    beam = HermiteBeam(w0, lambda, n, m);

                case 'laguerre'
                    beam = LaguerreBeam(w0, lambda, l, p);

                case {'elegant_hermite', 'eleganth', 'elegant_hg'}
                    beam = ElegantHermiteBeam(w0, lambda, n, m);

                case {'elegant_laguerre', 'elegantl', 'elegant_lg'}
                    beam = ElegantLaguerreBeam(w0, lambda, l, p);

                case {'hankel', 'hankel_laguerre'}
                    beam = HankelLaguerre(w0, lambda, l, p, htype);

                otherwise
                    error('BeamFactory:unknownType', ...
                        'Unknown beam type "%s". Supported: gaussian, hermite, laguerre, elegant_hermite, elegant_laguerre, hankel.', ...
                        type);
            end
        end

        function types = supportedTypes()
            % supportedTypes - Return cell array of all supported type names.
            types = {'gaussian', 'hermite', 'laguerre', ...
                     'elegant_hermite', 'elegant_laguerre', 'hankel'};
        end
    end

    % -----------------------------------------------------------------
    % Private helpers
    % -----------------------------------------------------------------
    methods (Static, Access = private)
        function val = getOpt(args, key, default)
            % getOpt - Extract value for 'key' from Name-Value list args.
            val = default;
            for i = 1:2:numel(args)-1
                if strcmpi(args{i}, key)
                    val = args{i+1};
                    return;
                end
            end
        end
    end
end
