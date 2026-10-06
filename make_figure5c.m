T_AMB = 300;
NF_dB = 3;
T_LNA = 10^(NF_dB/10)*T_AMB;
T_REC = 10^(NF_dB/10)*T_AMB;
G_TX_G_RX_dB = 40;
G_TX_G_RX = 10^(G_TX_G_RX_dB/10);
f_c_list = [800e6, 2.4e9, 5e9, 10e9];
B_list   = [1e6, 5e6, 20e6, 50e6];

% f_c_list = [433e6, 1e9 2.4e9, 5e9, 11e9];
% B_list   = [6e6, 10e6, 5e6, 20e6, 50e6];

% --- Match old figure style ---
fig = figure('Color','w', 'Position', [100, 100, 500, 450]); 
hold on; box on; grid on;
range = {};
R = {};
for i = 1:length(B_list)
    f_c = f_c_list(i);
    lam = 3e8/f_c;
    B_i = B_list(i);
    r_max = 3e8 / (2 * B_i);
    r_i = linspace(1, r_max, 2000);

    eta_i = 2 * sqrt(T_LNA*T_REC)/(T_AMB+T_REC) ...
            * G_TX_G_RX * lam^2 ./ (4*pi*r_i).^2 ...
            .* sinc(B_i*2*r_i/3e8);
    % eta_i = 2 * G_TX_G_RX * lam^2 ./ (4*pi*r_i).^2 .* sinc(B_i*2*r_i/3e8);

    R_i = eta_i.^2 .* B_i ./ (2*log(2));

    % Remove invalid (<=0) values for log plotting
    valid_idx = R_i > 0;
    R_i = R_i(valid_idx);
    r_i = r_i(valid_idx);
    range{i} = r_i;
    R{i} = R_i;

    loglog(r_i, R_i, 'LineWidth', 3.8, ...
        'DisplayName', ...
        "$f_c$ = " + num2str(f_c/1e9) + " GHz, B = " + num2str(B_i/1e6) + " MHz");
end

% --- Axes style to mirror old figure ---
set(gca, 'XScale', 'log', 'YScale', 'log', ...
         'FontSize', 14, 'LineWidth', 1.2, ...
         'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
         'TickDir', 'out');

xl = xlabel('Range (m)', 'FontSize', 20);
set(xl, 'Units','normalized');
xl.Position(2) = -0.065;   % adjust downward (try -0.10 to -0.15)
ylabel('$\overline{R}$ (b/s)', 'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'latex');

legend('Location', 'northeast', ...
       'FontSize', 18, ...
       'Box', 'off', ...
       'Interpreter', 'latex');
set(gca, 'FontSize', 20, 'LineWidth', 1.2, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
ylim([1, 1e6]);
xlim([1 600])
% title('Capacity vs Range for Different f_c and B', 'FontSize', 20, 'FontWeight', 'bold');
grid on; box on;
set(gca, 'FontName', 'Helvetica');
ax = gca;
ax.XScale = 'log';

ax.XTick = [1 100];          % 10^0 and 10^2
ax.XTickLabel = {'10^0','10^2'};
exportgraphics(gcf, 'out/Fig5c.png', 'Resolution', 300);

% save("Fig5c.mat", "range", "R");