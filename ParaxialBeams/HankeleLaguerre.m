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
