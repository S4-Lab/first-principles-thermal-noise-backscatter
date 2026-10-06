%High speed MATLAB codes
 
%This program use SCL with small list size as a "filter", i.e., under the
%same noise realization, If CRC-SCL with smaller L is correct, then CRC-SCL
%with larger L must be correct. You may not believe in this conjecture, but you can have
%a try, This is true with high probability (almost 1)
 
%This program has another accelerator. If in lower snr, CRC-SCL is correct, then at higher snr, 
%the same CRC-SCL must be correct.  You may not believe in this conjecture, but you can have
%a try, This is true with high probability (almost 1).
 
%I know above methods sound dangerous, but it is safe under following 4
%conditions
 
%1. The same code word
%2. The same AWGN noise realization with distribution N(0, 1) 
%3. CRC must be used (Only SCL is not permitted)
%4. The same code construction in all SNR range.
 
%This program satisfies above 4 conditions. You can verify the bler performance by comparisons with existing results. 
 
%Since MATLAB is not good at recursive function (too many parameters to be passed)
%I cancel the well-known recursiveCalcP() and recursiveCalcB() proposed by
%I. Tal
%Instead, 'For' function is used. You may argue that we can use objective
%oriented (OO) style. However, OO in matlab is also slow.
 
%Besides, the algorithms of following papers are provided.
 
%How to Construct Polar Codes
%Fast Successive-Cancellation Decoding of Polar Codes: Identification and Decoding of New Nodes
%beta-expansion: A Theoretical Framework for Fast and Recursive Construction of Polar Codes
 
%Above three algorithms are made by myself so the correctness is not guaranteed. 




clear
addpath('GA/')
% addpath('HowToConstructPolarCode/')
addpath('NodeProcess/')
% addpath('BECconstruction/')
% addpath('PolarizaedChannelsPartialOrder/')
%adding above folders will take round 2 seconds

design_epsilon = 0.32;
crc_length = 24;
[gen, det, g] = get_crc_objective(crc_length);
n = 12;
N = 2^n;
% K = round(N*0.085) + crc_length;
snr = 0.45;
coderate_vec = [0.2:-0.01:0.19];
list_vec = [64] ;  %row vec, you can write it like [1 4 16 32 ...]. The first element is always 1 for acceleration purpose. The ramaining elements are power of two.
max_runs = 5e3;
max_err = 100;
resolution = 10;%the results are shown per max_runs/resolution.

bler_vec = zeros(1, length(coderate_vec));
ber_vec = zeros(1, length(coderate_vec));
tic
parfor i = 1:length(coderate_vec)
    K = round(N*coderate_vec(i)) + crc_length;
    [bler, ber] = simulation_jida(N, K, design_epsilon, max_runs, max_err, resolution, snr, list_vec, gen, det, g, crc_length);
    bler_vec(i) = bler;
    ber_vec(i) = ber;
end
toc
plot(coderate_vec, ber_vec, '*')