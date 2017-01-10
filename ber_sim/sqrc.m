function h = sqrc(syms, beta)
% h = sqrc(syms, beta);
% Generate a Square-Root Raised Cosine Pulse
% 'syms' is 1/2 the length of srrc pulse in symbol durations
% 'beta' is the rolloff factor: beta=0 gives the sinc function
h = zeros(1, length(syms));
for n = 1 : length(syms)
    if syms(n) == 0
        h(n) = 1 + beta * (4 / pi - 1);
    elseif syms(n) == -1 / (4 * beta)
        h(n) = -beta * (2 / pi * cos(pi / (4 * beta) * (1+beta)) - cos(pi / (4 * beta) * (1 - beta)));
    elseif syms(n)==1/(4*beta)
        h(n) = -beta * (2 / pi * cos(pi / (4 * beta) * (1 + beta)) - cos(pi / (4 * beta) * (1 - beta)));
    else
        h(n) = (4 * beta * syms(n) * cos(pi * syms(n) * (1 + beta)) + sin(pi * syms(n) * (1-beta))) / ((1 - (4 * beta * syms(n)) ^ 2) * pi * syms(n));
    end
end
end