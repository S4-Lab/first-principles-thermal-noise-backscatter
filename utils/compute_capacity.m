function [cap,cap_ub] = compute_capacity(B, W, R)
S0 = 1;                 % noise PSD in off state (W/Hz);
S1 = (1+R)/(1-R)*S0;    % noise PSD in on state (W/Hz);
y0 = S0*W/B;
y1 = S1*W/B;
sigma02 = S0^2*W/B;
sigma12 = S1^2*W/B;
sigmamix2 = 1/2*sigma02 + 1/2*sigma12 + 1/4*(y0-y1)^2;
lim1 = y0 - 10*sqrt(sigma02);
lim2 = y1 + 10*sqrt(sigma12);
Hy = -quad(@(x) integrand(x, y0, sigma02, y1, sigma12), lim1, lim2);
Hyx = 1/(2*log(2)) + 1/4*log2(2*pi*sigma02) + 1/4*log2(2*pi*sigma12);
cap = (Hy - Hyx) * B;
cap_ub = R^2*W/2/log(2);
end

