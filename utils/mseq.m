function s = mseq(m, taps, seed)
% mseq  Generate an m-sequence (maximal-length PN) of length 2^m-1 in ±1.
%   s = mseq(m)                 uses a default primitive polynomial and seed
%   s = mseq(m, taps, seed)     taps is vector of tap positions (1=LSB)
%                               seed is 1×m binary, nonzero state
% Output s is 01 col vector.

    if nargin < 1, error('Provide m'); end
    % Default primitive tap sets for m=3..10 (1=LSB tap, m=MSB feedback)
    defaults = containers.Map('KeyType','double','ValueType','any');
    defaults(3) = [3 2];        % x^3 + x^2 + 1
    defaults(4) = [4 1];        % x^4 + x + 1
    defaults(5) = [5 2];        % x^5 + x^2 + 1
    defaults(6) = [6 1];        % x^6 + x + 1
    defaults(7) = [7 1];        % x^7 + x + 1
    defaults(8) = [8 6 5 4];    % x^8 + x^6 + x^5 + x^4 + 1
    defaults(9) = [9 4];        % x^9 + x^4 + 1
    defaults(10)= [10 3];       % x^10 + x^3 + 1
    defaults(11) = [11 2];             % x^11 + x^2 + 1
    defaults(12) = [12 6 4 1];         % x^12 + x^6 + x^4 + x + 1
    defaults(13) = [13 4 3 1];         % x^13 + x^4 + x^3 + x + 1
    defaults(14) = [14 5 3 1];         % x^14 + x^5 + x^3 + x + 1
    defaults(15) = [15 1];             % x^15 + x + 1
    defaults(16) = [16 5 3 2];         % x^16 + x^5 + x^3 + x^2 + 1
    defaults(17) = [17 3];             % x^17 + x^3 + 1
    defaults(18) = [18 7];             % x^18 + x^7 + 1
    defaults(19) = [19 5 2 1];         % x^19 + x^5 + x^2 + x + 1
    defaults(20) = [20 3];             % x^20 + x^3 + 1

    if nargin < 2 || isempty(taps), taps = defaults(m); end
    if nargin < 3 || isempty(seed), seed = ones(1,m); end
    if all(seed==0), error('Seed cannot be all zeros.'); end

    N = 2^m - 1;
    reg = seed(:).';  % shift register, reg(1)=LSB
    out = zeros(1,N);

    for k = 1:N
        out(k) = reg(end);                 % output MSB
        fb = 0;
        for t = 1:numel(taps)
            fb = xor(fb, reg(taps(t)));    % XOR taps (1-based LSB)
        end
        reg = [fb reg(1:end-1)];           % shift right, insert feedback at MSB
    end
    % s = 2*double(out)-1; % map {0,1} -> {-1,+1}
    s = out';
end