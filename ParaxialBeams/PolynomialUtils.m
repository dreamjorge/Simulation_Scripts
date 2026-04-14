classdef PolynomialUtils
    % PolynomialUtils - Special polynomial computations for paraxial beam modes
    % Compatible with GNU Octave and MATLAB
    %
    % Centralizes Hermite and Laguerre polynomial implementations so that
    % both standard and elegant beam classes share a single authoritative source.

    methods (Static)
        function H = hermitePoly(n, x)
            % hermitePoly - Physicist's Hermite polynomial H_n(x)
            % Recurrence: H_0=1, H_1=2x, H_{k+1}=2x*H_k - 2k*H_{k-1}
            if n == 0
                H = ones(size(x));
            elseif n == 1
                H = 2 * x;
            else
                H_prev2 = ones(size(x));
                H_prev1 = 2 * x;
                for k = 1:n-1
                    H = 2 * x .* H_prev1 - 2 * k * H_prev2;
                    H_prev2 = H_prev1;
                    H_prev1 = H;
                end
            end
        end

        function L = associatedLaguerre(p, l, x)
            % associatedLaguerre - Associated Laguerre polynomial L_p^|l|(x)
            % Explicit summation formula via binomial coefficients.
            n = p;
            k = abs(l);
            L = zeros(size(x));
            for m = 0:n
                L = L + ((-1)^m * nchoosek(n + k, n - m) .* x.^m) ./ factorial(m);
            end
        end
    end
end
