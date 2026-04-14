classdef HankelHermite
    % HankelHermite - Legacy-compatible Hankel-Hermite wrapper.
    %
    % Preserves historical constructor and static ray propagation APIs used
    % by research scripts:
    %   HankelHermite(x, y, hermiteParameters, hankelType)
    %   HankelHermite.getPropagateCartesianRays(...)

    properties
        HankelType
        OpticalField
    end

    methods (Static)
        function Rays = getPropagateCartesianRays(Rays, x, y, difr, HParametersZi, HParametersZ, HankelType)
            tempdr = num2cell(difr);
            [dx, dy, dz] = deal(tempdr{:});

            totalRays = numel(Rays.xCoordinate);

            for ray_index = 1:totalRays
                xi = Rays.xCoordinate(ray_index) + (1 ./ Rays.zxSlope(ray_index)) * dz;
                yi = Rays.yCoordinate(ray_index) + (1 ./ Rays.zySlope(ray_index)) * dz;
                zi = Rays.zCoordinate(ray_index) + dz;

                Rays.xCoordinate(ray_index) = xi;
                Rays.yCoordinate(ray_index) = yi;
                Rays.zCoordinate(ray_index) = zi;

                HHx = HankelHermite(x, yi, HParametersZi, HankelType);
                HHy = HankelHermite(xi, y, HParametersZi, HankelType);
                HHz = HankelHermite(xi, yi, HParametersZ, HankelType);

                fx = unwrap(angle(HHx.OpticalField));
                fy = unwrap(angle(HHy.OpticalField));
                fz = unwrap(angle(HHz.OpticalField));

                [zxSlope, zySlope, xySlope] = HankelHermite.gradientCartesian(fx, fy, fz, HParametersZi.k, dx, dy, dz, xi, yi, zi);
                Rays.zxSlope(ray_index) = zxSlope;
                Rays.zySlope(ray_index) = zySlope;
                Rays.xySlope(ray_index) = xySlope;
            end
        end
    end

    methods
        function obj = HankelHermite(x, y, hermiteParameters, hankelType)
            if nargin == 0
                obj.HankelType = 11;
                obj.OpticalField = [];
                return;
            end

            obj.HankelType = hankelType;
            waistGauss = hermiteParameters.Waist;

            [Hx, NHx] = HermiteParameters.getHermiteSolutions(hermiteParameters.n, (sqrt(2) ./ waistGauss) .* x);
            [Hy, NHy] = HermiteParameters.getHermiteSolutions(hermiteParameters.m, (sqrt(2) ./ waistGauss) .* y);

            gaussBeam = GaussianBeam(hermiteParameters.InitialWaist, hermiteParameters.Lambda);
            gaussField = gaussBeam.opticalField(x, y, hermiteParameters.zCoordinate);

            switch hankelType
                case 11
                    obj.OpticalField = (Hx + 1i * NHx) .* (Hy + 1i * NHy) .* gaussField;
                case 12
                    obj.OpticalField = (Hx + 1i * NHx) .* (Hy - 1i * NHy) .* gaussField;
                case 21
                    obj.OpticalField = (Hx - 1i * NHx) .* (Hy + 1i * NHy) .* gaussField;
                case 22
                    obj.OpticalField = (Hx - 1i * NHx) .* (Hy - 1i * NHy) .* gaussField;
                otherwise
                    error('HankelHermite:invalidType', 'Unsupported Hankel type: %d. Use 11, 12, 21, or 22.', hankelType);
            end
        end
    end

    methods (Static, Access = private)
        function [zxSlope, zySlope, xySlope] = gradientCartesian(fx, fy, fz, k, dx, dy, dz, x, y, z)
            gx = gradient(fx) ./ dx;
            gy = gradient(fy) ./ dy;
            gz = gradient(fz) ./ dz + k;

            n = size(gx, 2);
            idxZ = HankelHermite.clampIndex(floor(z / dz + 1), numel(gz));
            idxX = HankelHermite.clampIndex(n / 2 + 1 + floor(x / dx), numel(gx));
            idxY = HankelHermite.clampIndex(n / 2 + 1 + floor(y / dy), numel(gy));

            zxSlope = gz(idxZ) / gx(idxX);
            zySlope = gz(idxZ) / gy(idxY);
            xySlope = gx(idxX) / gy(idxY);
        end

        function idx = clampIndex(raw, maxSize)
            idx = floor(raw);
            if idx < 1
                idx = 1;
            elseif idx > maxSize
                idx = maxSize;
            end
        end
    end
end
