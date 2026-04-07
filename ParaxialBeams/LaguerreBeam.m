classdef LaguerreBeam
    % LaguerreBeam - Scalar optical field for a Laguerre-Gaussian beam
    % Compatible with GNU Octave and MATLAB
    
    properties
        Parameters      % LaguerreParameters object
        OpticalField    % 2D array representing the field
        r               % radial coordinate matrix
        theta           % angular coordinate matrix
    end
    
    methods
        function obj = LaguerreBeam(r, theta, params)
            % Constructor
            % r: radial coordinate matrix
            % theta: angular coordinate matrix
            % params: LaguerreParameters object
            
            if nargin > 0
                obj.Parameters = params;
                obj.r = r;
                obj.theta = theta;
                
                % Fundamental Gaussian Field
                GB = GaussianBeam(r, params);
                GField = GB.OpticalField;
                
                % Laguerre mode part
                w = params.Waist;
                l = params.l;
                p = params.p;
                
                % Normalization constant (optional but good for consistency)
                % sqrt(2*p! / (pi*(p+|l|)!))
                normConst = sqrt(2 * factorial(p) / (pi * factorial(p + abs(l))));
                
                % Amplitude term: (sqrt(2)*r/w)^|l|
                amplitudeTerm = (sqrt(2) * r ./ w).^abs(l);
                
                % Laguerre polynomial: L_p^|l|(2*r^2 / w^2)
                xArg = 2 * r.^2 ./ w.^2;
                Lpl = LaguerreParameters.getAssociatedLaguerrePolynomial(p, l, xArg);
                
                % Phase terms:
                % exp(i*l*theta)
                % exp(i*PhiPhase) -- where PhiPhase = (abs(l) + 2*p)*Gouy
                phaseMode = exp(1i * l * theta) .* exp(i * params.PhiPhase);
                
                obj.OpticalField = normConst .* (1./w) .* amplitudeTerm .* Lpl .* phaseMode .* GField;
                
                % Wait! GField already contains (w0/w). 
                % And GField contains exp(-r^2/w^2) and exp(i*k*z + i*k*r^2/2R - i*psi)
                % The standard LG formula usually collects everything.
                % If we multiply by GField, we must be careful not to double count terms.
                
                % Actually, GField = (w0/w) * exp(-r^2/w^2) * exp(-i*k*z + i*k*r^2/2R - i*psi)
                % The total phase shift for LG is (2p + |l| + 1)*psi.
                % Gaussian has 1*psi. So we add (2p + |l|)*psi. Correct.
                
                % One detail: LaguerreAmplitude in old code had (1/w).
                % Since GField has (w0/w), the total is (w0/w^2).
                % But the standard LG is usually normalized to total power 1.
                % I'll follow the old code's logic if possible, but keep it clean.
                
                % Updated calculation to avoid double counting:
                obj.OpticalField = amplitudeTerm .* Lpl .* exp(1i * l * theta) .* exp(1i * params.PhiPhase) .* GField;
            end
        end
    end
end