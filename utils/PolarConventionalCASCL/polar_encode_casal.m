function [x, frozen_bits, H_crc, lambda_offset, llr_layer_vec, bit_layer_vec, info_with_crc] = polar_encode_casal(N, K, info, crc_length)

lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

R = (K - crc_length)/N;
[gen, det, g] = get_crc_objective(crc_length);
[G_crc, H_crc] = crc_generator_matrix(g, K - crc_length);
crc_parity_check = G_crc(:, K - crc_length + 1 : end)';

design_snr = 2.5;
sigma_cc = 1/sqrt(2 * R) * 10^(-design_snr/20);
[channels, ~] = GA(sigma_cc, N);
[~, channel_ordered] = sort(channels, 'descend');
info_bits = sort(channel_ordered(1 : K), 'ascend');
frozen_bits = ones(N , 1);
frozen_bits(info_bits) = 0;
info_bits_logical = logical(mod(frozen_bits + 1, 2));

info_with_crc = [info; mod(crc_parity_check * info, 2)];
u = zeros(N, 1);
u(info_bits_logical) = info_with_crc;
x = polar_encoder(u, lambda_offset, llr_layer_vec);
end

