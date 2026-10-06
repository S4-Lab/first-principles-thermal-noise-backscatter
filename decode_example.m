%% TX
addpath utils/
addpath utils/PolarConventionalCASCL/;
addpath utils/PolarConventionalCASCL/GA/;
addpath utils/PolarConventionalCASCL/NodeProcess/;

data = load("rawmeasurements/example_decode_data.mat");

bw = 50e6;
fb = 10000;
preamble = mseq(11);
preamble_bits_num = length(preamble);
K = data.K;
E = data.E;
msg = data.msg;
payload_bits_num = E;
bit_samples = round(bw / fb);
crc_len = 24;
[payload, frozen_bits, H_crc, lambda_offset, llr_layer_vec, bit_layer_vec, msg_crc] = polar_encode_casal(E, K, msg, crc_len);

%% RX - Synchronization
% received signal from USRP
powers = data.powers';
% sync
filter_kernal = ones(preamble_bits_num * bit_samples, 1);
for i = 1:preamble_bits_num
    if preamble(i) == 0
        filter_kernal((i-1)*bit_samples+1:i*bit_samples) = ...
            -1*filter_kernal((i-1)*bit_samples+1:i*bit_samples);
    end
end
dsamp = 1;
[c, l] = xcorr(powers, filter_kernal);
valid_idx = l >= 0 & l <= length(powers) - length(filter_kernal);
c = c(valid_idx);
c = c - mean(c);
l = l(valid_idx);
[~, start_idx] = max(abs(c));
% check bit flip
bit_flip = c(start_idx) < 0;
% plot
plot_dsamp = 1000;
figure; hold on;
plot(l(1:plot_dsamp:end) / bw, c(1:plot_dsamp:end), LineWidth = 2.5);
xline(start_idx / dsamp / bw, 'r--', LineWidth=1.5);
xlabel('Time (sec)');
ylabel('Correlation')

%% RX - Channel estimation
powers_preamble = powers(start_idx + 1 : start_idx + preamble_bits_num * bit_samples);
preamble_demodulate = mean(reshape(powers_preamble, bit_samples, preamble_bits_num), 1)';
preamble_demodulate_1 = preamble_demodulate(preamble == 1);
preamble_demodulate_0 = preamble_demodulate(preamble == 0);
T1 = mean(preamble_demodulate_1);
T0 = mean(preamble_demodulate_0);
sigma2 = (var(preamble_demodulate_1) + var(preamble_demodulate_1)) / 2;
threshold = (T1 + T0) / 2;
BER_emp = (sum(preamble_demodulate_1 < threshold)...
    + sum(preamble_demodulate_0 > threshold))/ preamble_bits_num;
if bit_flip
    BER_emp = 1 - BER_emp;
end
R_emp = (T1 - T0) / (T1 + T0);
fprintf("Empirical R: %.5f\n", R_emp);
fprintf("Empirical uncoded BER: %.5f\n", BER_emp);
figure; hold on;
histogram(preamble_demodulate_0, DisplayName="0-bit");
histogram(preamble_demodulate_1, DisplayName="1-bit");
xline(threshold, 'b--', LineWidth=3, DisplayName="Thresh")
xlabel("Received Power"); ylabel("Count"); legend(FontSize=24);
%% RX - Payload Decode
payload_idx = start_idx + preamble_bits_num*bit_samples;
powers_payload = powers(payload_idx + 1 : payload_idx + payload_bits_num * bit_samples);
payload_demodulate = mean(reshape(powers_payload, bit_samples, payload_bits_num), 1)';
LLR = 1/2/sigma2 * ((payload_demodulate - T1).^2 - (payload_demodulate - T0).^2); % log(pr(0)/pr(1))
fprintf("Payload BER Before ECC: %.5f\n", sum((LLR < 0 )~=payload) / payload_bits_num)

list_size = 1024;
decBits = CASCL_decoder(LLR, list_size, K, frozen_bits, H_crc, lambda_offset, llr_layer_vec, bit_layer_vec);
fprintf("Decoded Info Bits BER: %.5f\n", 1-sum(decBits == msg_crc) / K);
fprintf('Decoded Bits: %s...\n', sprintf('%d', decBits(1:10)));
