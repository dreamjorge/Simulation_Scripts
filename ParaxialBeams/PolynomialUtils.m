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

        function Xnmr = xAssociatedLaguerre(n, m, x)
            % xAssociatedLaguerre - Second independent solution of the associated
            % Laguerre differential equation.
            %
            % Three-term series expansion involving logarithmic and digamma (psi)
            % functions. Ported from dreamjorge/LaguerreGaussBeams.
            %
            % References:
            %   - Papi, J. "Hankel beams", Structured Light, Elsevier 2023.

            n_terms = 97;

            % Term 1: power series for k = n+2 .. n_terms
            a1 = x.^(n+1) ./ (gamma(n+2) * gamma(m+1+n+1));
            for k = (n+2):n_terms
                a1 = a1 + gamma(k-n) .* x.^k ./ (gamma(k+1) * gamma(m+1+k));
            end

            % Term 2: negative power series for k = 2..m (only when m >= 1)
            a2 = zeros(size(x));
            if m >= 1
                a2 = x.^(-1) ./ (gamma(n+2) * gamma(m));
                for k = 2:m
                    a2 = a2 + gamma(k) .* x.^(-k) ./ (gamma(n+k+1) * gamma(m-k+1));
                end
            end

            % Term 3: logarithmic + digamma series
            a3 = (-log(x) + psi(1) + psi(m+1) - psi(n+1)) ...
                 ./ (gamma(n+1) * gamma(m+1));
            for k = 1:n
                a3 = a3 + (-log(x) + psi(k+1) + psi(m+k+1) - psi(n-k+1)) ...
                     .* (-x).^k ./ (gamma(n-k+1) * gamma(m+k+1) * gamma(k+1));
            end

            Xnmr = (gamma(n+m+1) / pi) .* (-a1 + (-1)^n .* (a2 + a3));
        end

        function Xn = xLaguerreG(n, m, x)
            % xLaguerreG - Full radial envelope for XLG (second Laguerre solution).
            %
            % Combines the second associated Laguerre solution with the
            % standard exponential-polynomial envelope.
            %
            % IMPORTANT:
            %   This returns a FULL envelope term:
            %     exp(-x/2) * x^(m/2) * xAssociatedLaguerre(...)
            %   Use xAssociatedLaguerre() directly when beam assembly already
            %   applies Gaussian carrier and (sqrt(2)r/w)^m scaling.
            %
            % Ported from dreamjorge/LaguerreGaussBeams.

            Xn = (-1)^(n+1) ./ ((n + (m+1)/2).^(m/2)) ...
                 .* exp(-x./2) .* x.^(m/2) ...
                 .* PolynomialUtils.xAssociatedLaguerre(n, m, x);
            Xn(isnan(Xn) | isinf(Xn)) = 0;
        end
    end
end
