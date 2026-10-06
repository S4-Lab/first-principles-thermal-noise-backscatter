% Make sure to set the current folder to be ThermalNoiseCode
%% plot capacity result
bw = 50e6;
RR = (0.001:0.00001:0.015)';
dist = [];

fig = figure('Color','w', 'Position', [100, 100, 500, 450]);
hold on; box on; grid on;

% --- Theoretical curve ---
C_th = RR.^2 * bw / (2 * log(2));
loglog(RR, C_th, 'k-', 'LineWidth', 3.5, 'DisplayName', 'Bound $\overline{R}$');

% --- Read measured data ---
rootDir = 'rawmeasurements/R_casal/';
files = dir(fullfile(rootDir, '**', '*'));
files = files(~[files.isdir]);
files = files(~startsWith({files.name}, '.'));

R_list = [];
dtrate_list = [];
fb_list = [];
gap_list = [];
finite_fb_gap_list = [];

for i = 1:length(files)
    pp = fullfile(files(i).folder, files(i).name);
    parts = split(string(pp), '/');
    fb = str2double(parts{end-1});
    R_i = load(pp, "R_emp");
    E_i = load(pp, "E");
    K_i = load(pp, "K");
    dr = K_i.K / E_i.E * fb;
    R_list = [R_list ; R_i.R_emp];
    dtrate_list = [dtrate_list ; dr];
    fb_list = [fb_list ; fb];
    dist = [dist ; dr / fb];

    snr_emp = 2^(2*dr/fb)-1;
    snr = bw/fb*(R_i.R_emp)^2;
    Gamma = snr_emp / snr;
    [cap_theory, ~] = compute_capacity(fb, bw, R_i.R_emp);

    gap = dr / ((R_i.R_emp)^2*bw/2/log(2));
    
    gap_theory = cap_theory / ((R_i.R_emp)^2*bw/2/log(2));
    % gap_list = [gap_list ; gap];
    gap_list = [gap_list ; Gamma];
    finite_fb_gap_list = [finite_fb_gap_list ; gap_theory];
end

% --- Unique fb values and color mapping ---
fb_values = [500, 1000, 2500, 5000, 10000, 20000, 40000];
colors = jet(length(fb_values));  % jet goes from blue (cold) to red (hot)

% --- Scatter data by fb ---
for k = 1:length(fb_values)
    idx = fb_list == fb_values(k);
    if any(idx)
        scatter(R_list(idx), dtrate_list(idx), 70, ...
            'filled', ...
            'MarkerFaceColor', colors(k,:), ...
            'MarkerEdgeColor', 'k', ...
            'DisplayName', sprintf('$R_{b} = %.1f$ kb/s', fb_values(k)/1000));
    end
end

% --- Axes setup ---
set(gca, 'XScale', 'log', 'YScale', 'log');
xl = xlabel('Power Contrast Ratio $\widetilde{\Delta P}$', ...
    'FontSize', 28, 'FontWeight', 'bold', 'Interpreter', 'latex');
set(xl, 'Units','normalized');
xl.Position(2) = -0.07;   % adjust downward (try -0.10 to -0.15)
ylabel('$R_{\bf emp}$ (b/s)', 'FontSize', 28, 'FontWeight', 'bold', 'Interpreter', 'latex');
set(gca, 'FontSize', 20, 'LineWidth', 1.2, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
legend('Location', 'southeast', 'FontSize', 20, 'Box', 'off', Interpreter='latex');
ylim([10 1e4])
xlim([0.001 0.013])
% title('Empirical vs. Theoretical Channel Capacity', ...
    % 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontName', 'Helvetica');
exportgraphics(gcf, 'out/Fig5a.png', 'Resolution', 300);
% power_contrast_ratio_bound = RR;
% R_bound = C_th;
% power_contrast_ratio_empirical = R_list;
% R_empirical = dtrate_list;
% save("Fig5a.mat", "power_contrast_ratio_bound", "R_bound", "power_contrast_ratio_empirical", "R_empirical");
%% Plot capacity vs range (log scale)
fig = figure('Color','w', 'Position', [100, 100, 500, 450]); hold on; box on; grid on;

% --- Define root directories ---
rootDir = { ...
    'rawmeasurements/range/'
};

% --- Loop through each case ---
for jj = 1:numel(rootDir)
    rootd = rootDir{jj};

    % Get all data files (ignore folders and hidden files)
    files = dir(fullfile(rootd, '**', '*'));
    files = files(~[files.isdir]);
    files = files(~startsWith({files.name}, '.'));

    ranges_ = [];
    capas_  = [];
    Rs_ = [];

    % --- Read each file ---
    for i = 1:numel(files)
        pp = fullfile(files(i).folder, files(i).name);
        parts = split(string(pp), filesep);   % use filesep for cross-platform paths

        % extract range (assuming it's always in part 10)
        r = str2double(parts{end-1});

        if isnan(r)
            continue
        end

        % load datarate from .mat file
        data = load(pp, "datarate");
        if isfield(data, "datarate")
            cc = data.datarate;
        else
            continue;
        end
        data = load(pp, "R_emp");
        R_emp_ = data.R_emp;

        ranges_ = [ranges_; r];
        capas_  = [capas_; cc];
        Rs_ = [Rs_ ; R_emp_];
    end

    % --- Clean and sort data ---
    valid = capas_ > 0 & ~isnan(capas_);
    ranges_ = ranges_(valid);
    capas_  = capas_(valid);

    [ranges_, idx] = sort(ranges_);
    capas_ = capas_(idx);

    % --- Plot using semilogy ---
    semilogy(ranges_+2, capas_, '^-', 'LineWidth', 2.8, DisplayName="Empirically Achieved");
    semilogy(ranges_+2, Rs_.^2*bw/2/log(2), 'black--', 'LineWidth', 2.8, DisplayName="Bound $\overline{R}$");
end

% --- Final formatting ---
xl = xlabel('Range (cm)', 'FontSize', 16);
set(xl, 'Units','normalized');
xl.Position(2) = -0.07;   % adjust downward (try -0.10 to -0.15)
ylabel('$R_{\rm emp}$ (b/s)', 'FontSize', 16, 'Interpreter', 'latex');

set(gca, 'YScale', 'log'); % ensure semilogy
set(gca, 'FontSize', 20, 'LineWidth', 1.2, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
legend('Location', 'northeast', 'FontSize', 20, 'Box', 'off', Interpreter='latex');
grid on; box on;
xlim([25 90])
% title('Measured Capacity vs Range', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontName', 'Helvetica');
exportgraphics(gcf, 'out/Fig5b.png', 'Resolution', 300);

% range = ranges_;
% R_bound = Rs_.^2*bw/2/log(2);
% R_empirical = capas_;
% save("Fig5b.mat", "range", "R_bound", "R_empirical");