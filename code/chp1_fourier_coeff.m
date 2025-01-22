% Calculates Fourier Coefficients using numerical integration

function [a,b] = FourierCoeff(f,T,n)
    % Compute Fourier coefficients a_n and b_n for given function f
    % Define integration function handles
    cos_term = @(x) f(x) .* cos(2 * pi * n * x / T);
    sin_term = @(x) f(x) .* sin(2 * pi * n * x / T);

    % Perform numerical integration
    a = integral(cos_term, 0, T) * 2 / T;
    b = integral(sin_term, 0, T) * 2 / T;
end