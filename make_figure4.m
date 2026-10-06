clear all; clc;
addpath utils/;
% system params
bw = 50e6;
function [powers_, bits_0_idx, bits_1_idx, bits] = sys_sim(fb, R)
    
    rng(43, 'twister');

    T1 = 1+R;
    T0 = (1-R) / (1+R) * T1;
    bw = 50e6;
    q = sqrt(bw / fb) * R;
    BER = qfunc(q)
    [capacity_achieve, capacity] = compute_capacity(fb, bw, R)
    
    % simulate 100000 bits
    bits = randi([0, 1], [1000000, 1]);
    bits(1:256) = repmat([0;0;1;1], [64, 1]);
    
    bits_0_idx = bits == 0;
    bits_1_idx = bits == 1;
    sigma02 = T0^2*fb/bw;
    sigma12 = T1^2*fb/bw;
    
    powers_ = zeros(size(bits));
    powers_(bits_0_idx) = T0 + sqrt(sigma02) * randn(sum(bits_0_idx), 1);
    powers_(bits_1_idx) = T1 + sqrt(sigma12) * randn(sum(bits_1_idx), 1);
end

%% Plot time domain
% --- Configuration & Parameters ---
fb_values = [1000, 5000, 25000];
R = 2e-2;
bits_2_plot = 20;

% Total height 1200 (3 plots x 400) to keep your original aspect ratio
% figure('Color', 'w', 'Position', [100, 100, 700, 1200]);
figure('Color', 'w', 'Position', [100, 100, 2400, 450]);

% Color Definitions & Lightening
c_blue = [0, 0.4470, 0.7410]; % Demodulated Blue
c_red  = [0.85, 0.3, 0.1];    % Binary Red
lighten = @(color, factor) color + (1 - color) * factor;

for k = 1:length(fb_values)
    subplot(1, 3, k);
    hold on; box on;
    
    fb = fb_values(k);
    
    % --- Constants & Simulation ---
    T1 = 1+R;
    T0 = (1-R) / (1+R) * T1;
    [powers_, bits_0_idx, bits_1_idx] = sys_sim(fb, R); % Ensure sys_sim is in path
    
    % Time Vector
    t_plot = (0:bits_2_plot-1)'/fb*1000; 
    
    % 1. Binary Theoretical (Shifted Square Wave)
    p_theoretical = T0 * ones(bits_2_plot, 1);
    p_theoretical(bits_1_idx(1:bits_2_plot)) = T1;
    
    t_stairs = ((0:bits_2_plot)' - 0.5) / fb * 1000; 
    p_stairs = [p_theoretical; p_theoretical(end)]; 
    
    % Using lightened red for consistency with previous figure style
    stairs(t_stairs, p_stairs, '--', ...
        'Color', lighten(c_red, 0.5), ...  
        'LineWidth', 4, ... % Reduced slightly for stack readability
        'DisplayName', 'Binary Bits');
        
    % 2. Demodulated
    powers_plot = powers_(1:bits_2_plot);
    plot(t_plot, powers_plot, ...
        'Color', lighten(c_blue, 0.3), ... % Lightened blue
        'Marker', 'o', ...
        'MarkerSize', 8, ...
        'LineWidth', 3, ... 
        'DisplayName', 'Demodulated Y');

    % --- Nature Style Formatting ---
    set(gca, ...
        'FontName', 'Helvetica', ... 
        'FontSize', 28, ... 
        'LineWidth', 1.2, ...      
        'TickDir', 'out');       
    
    % Dynamic X-limit based on fb
    xlim([0, (bits_2_plot-1.5)/fb*1000]);
    ylim([0.94, 1.06]); % Tightened to show signal variations clearly
    
    xlabel('t (ms)', 'FontName', 'Helvetica', 'FontSize', 28);
    ylabel('Norm. Power', 'FontName', 'Helvetica', 'FontSize', 28);
    
    % Legend only on the first plot to save space, or all if preferred
    if k == 1
        legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 24);
    end
    
    hold off;
    % T = t_plot;
    % P = powers_plot;
    % save("Fig4b_"+num2str(k)+".mat", "T", "P")
end

% Save the final figure
exportgraphics(gcf, 'out/Fig4b.png', 'Resolution', 300);

%% Plot distribution
% --- Configuration & Parameters ---
fb_values = [1000, 5000, 25000];
R = 2e-2;

% Figure height 1200 (3 plots x 400) to match your previous layouts
% figure('Color', 'w', 'Position', [950, 100, 700, 1200]);
figure('Color', 'w', 'Position', [100, 100, 2400, 450]);

% Color Definitions (Matching your previous lighter blue and orange/red)
c_bit0 = [0.1, 0.3, 0.7];   % Blue 
c_bit1 = [0.85, 0.3, 0.1];  % Red/Orange
lighten = @(color, factor) color + (1 - color) * factor;

for k = 1:length(fb_values)
    subplot(1, 3, k);
    hold on; box on;
    
    fb = fb_values(k);
    
    % --- Simulation Block ---
    [powers_, bits_0_idx, bits_1_idx, bits] = sys_sim(fb, R);
    
    % Separate powers based on bit value
    powers_0 = powers_(bits == 0);
    powers_1 = powers_(bits == 1);
    
    % --- Plot Histograms ---
    % Using FaceAlpha 0.5 as requested for transparency
    histogram(powers_0, 50, ...
        'Normalization', 'probability', ...
        'FaceColor', lighten(c_bit0, 0.3), ... 
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.5, ...               
        'DisplayName', 'Bit 0');

    histogram(powers_1, 50, ...
        'Normalization', 'probability', ...
        'FaceColor', lighten(c_bit1, 0.3), ...  
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.5, ...               
        'DisplayName', 'Bit 1');

    % --- Threshold Logic ---
    threshold = (mean(powers_0) + mean(powers_1)) / 2; 
    xline(threshold, '--k', ...
        'LineWidth', 2, ...
        'DisplayName', 'Threshold');

    % --- Nature Style Formatting ---
    set(gca, ...
        'FontName', 'Helvetica', ... 
        'FontSize', 28, ... 
        'LineWidth', 1.2, ...      
        'TickDir', 'out');       

    xlim([0.9 1.1]); % Tightened to focus on the overlap area
    ylim([0 0.1]);     % Adjusted slightly for vertical stacking
    
    xlabel('Norm. Power', 'FontName', 'Helvetica', 'FontSize', 28);
    ylabel('PDF', 'FontName', 'Helvetica', 'FontSize', 28);

    % Legend only on the first plot to maximize space
    if k == 1
        legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 24);
    end
    
    hold off;
    % P_0 = powers_0; P_1 = powers_1;
    % save("Fig4c_"+num2str(k)+".mat", "P_0", "P_1")
end

% Save the final figure
exportgraphics(gcf, 'out/Fig4c.png', 'Resolution', 300);


%% Plot BER and Capacity

% Define vectors and definitions that were missing in snippet
fb_vec = linspace(100, 200000, 1000)'; % Using User's range (Linear)

% Re-calculate BER for this specific fb_vec
q_vec_sweep = sqrt(bw ./ fb_vec) * R;
BER_vec = 0.5 * erfc(q_vec_sweep/sqrt(2));

% Define markers and colors (Reduced range to match 100-10000 Hz)
marker_fbs = [1000, 5000, 25000]; 
marker_colors = {[0.47, 0.67, 0.19], [0.85, 0.3, 0.1], [0.5, 0.1, 0.5]};

capacity_achieve_vec = zeros(length(fb_vec), 1);
for i = 1:length(fb_vec)
    [capacity_achieve_vec(i), capacity_limit] = compute_capacity(fb_vec(i), bw, R);
end

limit_kbps = capacity_limit / 1000;

figure('Color', 'w', 'Position', [750 1200 1600 350]); % Increased height slightly for larger fonts
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---------- Tile 1: Achieved Rate ----------
ax1 = nexttile(2); hold(ax1, 'on'); box on;
semilogx(fb_vec, capacity_achieve_vec / 1000, 'Color', [0.1 0.3 0.7], 'LineWidth', 3);
% semilogx(fb_vec, 1/2*fb_vec.*log2(1+bw*R^2./fb_vec) / 1000, 'Color', [1 0 0], 'LineWidth', 3);
yline(capacity_limit / 1000, 'Color', [0.1 0.1 0.1], ...
    'LineWidth', 5, ...
    'Label', 'Capacity Upper Bound', ... % Use \bar for tex
    'FontSize', 20, ...
    'FontName', 'Helvetica', ...
    'LabelVerticalAlignment', 'top', ...
    'LabelHorizontalAlignment', 'right', ...
    'Interpreter', 'tex'); % Change to tex
yline(capacity_limit / 1000, 'Color', [0.1 0.1 0.1], ...
    'LineWidth', 5, ...
    'Label', '14.4 kb/s', ... % Use \bar for tex
    'FontSize', 18, ...
    'FontName', 'Helvetica', ...
    'LabelVerticalAlignment', 'top', ...
    'LabelHorizontalAlignment', 'left', ...
    'Interpreter', 'tex'); % Change to tex

for k = 1:length(marker_fbs)
    target_f = marker_fbs(k);
    [~, idx] = min(abs(fb_vec - target_f));
    val = capacity_achieve_vec(idx) / 1000;
    
    % Vertical and Horizontal Intercepts
    xline(target_f, ':', 'Color', marker_colors{k}, 'LineWidth', 2);
    yline(val, '--', 'Color', marker_colors{k}, 'LineWidth', 1.5, ...
        'Label', sprintf('%.1f kb/s', val), 'FontSize', 18, ... % Increased Intercept Font
        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');

     % 2. Draw Double-Sided Arrow (from dot to capacity bound)
    % We use 'quiver' to create arrowheads at both ends
    arrow_head_size = 3;
    quiver(target_f, val, 0, limit_kbps - val, 0, ...
        'Color', marker_colors{k}, 'LineWidth', 2, 'MaxHeadSize', arrow_head_size);
    quiver(target_f, limit_kbps, 0, val - limit_kbps, 0, ...
        'Color', marker_colors{k}, 'LineWidth', 2, 'MaxHeadSize', arrow_head_size);

    % 3. Calculate and Plot Gap Text
    gap_val = limit_kbps - val;
    gap_percent = (1 - (val/limit_kbps)) * 100;
    
    % Placing text to the right of the arrow
    text(target_f * 1.1, limit_kbps*0.93, ...
        sprintf('Gap: %.1f%%', gap_percent), ...
        'Color', marker_colors{k}, 'FontSize', 18, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left');
    
    
    plot(target_f, val, 'o', 'MarkerSize', 10, 'MarkerFaceColor', marker_colors{k});
    
    % Offset text further away from the marker (using 1.15 multiplier)
    text(target_f*0.65, val + 1, sprintf('%d kb/s', target_f/1000), ...
        'Color', marker_colors{k}, 'FontSize', 22, ... % Increased Marker Font
        'HorizontalAlignment', 'center');
end
set(ax1, 'FontName', 'Arial', 'FontSize', 20, 'LineWidth', 1.2, 'XScale', 'log');
ylabel('R (kb/s)', 'Interpreter', 'tex', 'FontSize', 22);
xlabel('R_b (b/s)', 'Interpreter', 'tex', 'FontSize', 22);
ylim([0 max(capacity_achieve_vec/1000)*1.25]); % Increased Y-limit to make room for text
xlim([fb_vec(1) fb_vec(end)])

% ---------- Tile 2: BER ----------
ax2 = nexttile(1); hold(ax2, 'on'); box on;
semilogx(fb_vec, BER_vec, 'Color', [0.1 0.3 0.7], 'LineWidth', 3);

for k = 1:length(marker_fbs)
    target_f = marker_fbs(k);
    [~, idx] = min(abs(fb_vec - target_f));
    val_ber = BER_vec(idx);
    
    % Vertical and Horizontal Intercepts
    xline(target_f, ':', 'Color', marker_colors{k}, 'LineWidth', 2);
    yline(val_ber, '--', 'Color', marker_colors{k}, 'LineWidth', 1.5, ...
        'Label', sprintf('%.2f', val_ber), 'FontSize', 18, ... % Increased Intercept Font
        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');
    
    plot(target_f, val_ber, 'o', 'MarkerSize', 10, 'MarkerFaceColor', marker_colors{k});
    
    % Offset text further away from the marker (fixed 0.04 offset)
    text(target_f*0.65, val_ber + 0.03, sprintf('%d kb/s', target_f/1000), ...
        'Color', marker_colors{k}, 'FontSize', 22, ... % Increased Marker Font
        'HorizontalAlignment', 'center');
end
set(ax2, 'FontName', 'Arial', 'FontSize', 20, 'LineWidth', 1.2, 'XScale', 'log');
ylabel('BER', 'FontSize', 22);
xlabel('R_b (b/s)', 'Interpreter', 'tex', 'FontSize', 22);
ylim([0 0.4]); % Updated to 0.4 as requested
xlim([fb_vec(1) fb_vec(end)])
set(gca, 'FontName', 'Helvetica');
linkaxes([ax1 ax2], 'x');
exportgraphics(gcf, 'out/Fig4de.png', 'Resolution', 300);

% Rb = fb_vec;
% BER = BER_vec;
% R = capacity_achieve_vec;
% save("Fig4d", "Rb", "BER");
% save("Fig4e", "Rb", "R");
