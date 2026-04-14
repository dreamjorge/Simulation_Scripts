classdef HankeleHermite
    % HankeleHermite - Legacy alias for HankelHermite.

    properties
        HankelType
        OpticalField
    end

    methods (Static)
        function Rays = getPropagateCartesianRays(Rays, x, y, difr, HParametersZi, HParametersZ, HankelType)
            Rays = HankelHermite.getPropagateCartesianRays(Rays, x, y, difr, HParametersZi, HParametersZ, HankelType);
        end
    end

    methods
        function obj = HankeleHermite(x, y, hermiteParameters, hankelType)
            % Deprecated: Use HankelHermite instead.
            % This alias will be removed in a future release.
            persistent warnIssued
            if isempty(warnIssued)
                warning('HankeleHermite is deprecated. Use HankelHermite instead. See docs/migration/LEGACY_MIGRATION_PLAN.md');
                warnIssued = true;
            end

            if nargin == 0
                obj.HankelType = 11;
                obj.OpticalField = [];
                return;
            end

            hankel = HankelHermite(x, y, hermiteParameters, hankelType);
            obj.HankelType = hankel.HankelType;
            obj.OpticalField = hankel.OpticalField;
        end
    end
end
