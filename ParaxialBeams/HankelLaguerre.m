classdef HankelLaguerre
   % HankelLaguerre - Hankel-type Laguerre-Gaussian beam field
    % Used in Hankel-based ray tracing (HankelLaguerrePropagation.m, etc.)
    %
    % Theoretical Context (Contexto Academico):
    % Hankel beams are specialized solutions to the paraxial wave equation 
    % that represent "conical" waves. Unlike standard LG modes, these beams 
    % possess a net radial energy flow, making them ideal for modeling 
    % diverging or converging wavefronts in complex media.
    %
    % Mathematical Definition:
    % They are constructed through the superposition of a paraxial mode 
    % and its quadrature (Hilbert transform) companion:
    %   $$ H_{lp}^{(1,2)} (r, \phi, z) = LG_{lp} (r, \phi, z) \pm i \cdot XLG_{lp} (r, \phi, z) $$
    %
    % Where XLG is the Hilbert-transformed counterpart with:
    %   - Phase: exp(-1i*p*theta) instead of exp(1i*l*theta)
    %   - Polynomial: AssociatedLaguerre(p, l, xArg) (same polynomial, different phase)
    %
    % Physical Significance:
    % - $H^{(1)}$: Outward propagating wave (away from axis).
    % - $H^{(2)}$: Inward propagating wave (towards axis).

    properties
        OpticalFieldLaguerre  % Complex field array
    end

    methods
        function obj = HankelLaguerre(r, theta, params, hankelType)
            % Constructor
            % r: radial coordinate matrix
            % theta: angular coordinate matrix
            % params: LaguerreParameters object
            % hankelType: 1 (H^(1)) or 2 (H^(2))
            
            if nargin < 4
                hankelType = 1;
            end
            
            obj.OpticalFieldLaguerre = computeHankelField(r, theta, params, hankelType);
        end
    end
end

function field = computeHankelField(r, theta, params, hankelType)
    l = params.l;
    p = params.p;
    w = params.Waist;
    
    % Standard LG amplitude term: (sqrt(2)*r/w)^|l|
    amp = (sqrt(2) * r ./ w).^abs(l);
    
    % x argument for Laguerre polynomial
    xArg = 2 * r.^2 ./ w.^2;
    
    % Standard LG polynomial L_p^l(x)
    Lpl = PolynomialUtils.associatedLaguerre(p, l, xArg);
    
    % Gaussian carrier field
    GB = GaussianBeam(r, params);
    GField = GB.OpticalField;
    
    % Standard LG field (LB)
    LB_field = amp .* Lpl .* exp(1i * l * theta) .* exp(1i * params.PhiPhase) .* GField;
    
    % XLG (Hilbert-transformed) field:
    % - Phase uses exp(-1i*p*theta) instead of exp(1i*l*theta)
    % - Polynomial is the same L_p^l(x)
    XLG_field = amp .* Lpl .* exp(-1i * p * theta) .* exp(1i * params.PhiPhase) .* GField;
    
    % Combine via Hankel type
    if hankelType == 1
        field = LB_field + 1i * XLG_field;
    else
        field = LB_field - 1i * XLG_field;
    end
end