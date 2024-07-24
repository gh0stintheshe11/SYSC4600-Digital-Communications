function h = root_raised_cosine(beta, L, T, eta)
%ROOT_RAISED_COSINE Root raised cosine filter impulse response
%
%   h = root_raised_cosine(beta, L, T, eta) gives the impulse response
%   of a root raised cosine filter with a roll-off factor of beta and
%   a symbol period of T, truncated to L symbols, with an oversampling
%   ratio of eta samples per symbol.

    t = (-L/2*eta:L/2*eta-1) / eta;

    bt = 4 * beta * t;

    h1 = bt .* cos(pi*t*(1+beta)) + sin(pi*t*(1-beta));
    h2 = pi*t .* (1 - bt.*bt);
    h = h1 ./ h2;

    h(t == 0) = 1 - beta + 4*beta/pi;
    h(abs(bt) == 1) = beta/sqrt(2) * ((1+2/pi)*sin(pi/(4*beta)) ...
                            + (1-2/pi)*cos(pi/(4*beta)));

    h = h ./ sqrt(T);

    h = h ./ sqrt(sum(h.*h)*T/eta);     % Normalize truncated pulse
end
