function answer = integrand(z, y0, sigma02, y1, sigma12)

py0 = 1/sqrt(2*pi*sigma02)*exp(-((z-y0).^2)/(2*sigma02));
py1 = 1/sqrt(2*pi*sigma12)*exp(-((z-y1).^2)/(2*sigma12));
py = 1/2*(py0 + py1);
answer = py.*log2(py);
end