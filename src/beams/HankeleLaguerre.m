classdef HankeleLaguerre
    % HankeleLaguerre - Legacy alias for HankelLaguerre.

    properties
        HankelType
        OpticalFieldLaguerre
    end

    methods (Static)
        function Rays = getPropagateCylindricalRays(Rays, TotalRays, r, th, difr, LParametersZi, LParametersZ, HankelType)
            Rays = HankelLaguerre.getPropagateCylindricalRays(Rays, TotalRays, r, th, difr, LParametersZi, LParametersZ, HankelType);
        end
    end

    methods
        function obj = HankeleLaguerre(rCoordinate, thetaCoordinate, laguerreParameters, hankelType)
            % Deprecated: Use HankelLaguerre instead.
            % This alias will be removed in a future release.
            persistent warnIssued
            if isempty(warnIssued)
                warning('HankeleLaguerre is deprecated. Use HankelLaguerre instead. See docs/migration/LEGACY_MIGRATION_PLAN.md');
                warnIssued = true;
            end

            if nargin == 0
                obj.HankelType = 1;
                obj.OpticalFieldLaguerre = [];
                return;
            end

            hankel = HankelLaguerre(rCoordinate, thetaCoordinate, laguerreParameters, hankelType);
            obj.HankelType = hankel.HankelType;
            obj.OpticalFieldLaguerre = hankel.OpticalFieldLaguerre;
        end
    end
end
