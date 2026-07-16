function analyze_snow_dependence(repoRoot)
%
% ANALYZE_SNOW_DEPENDENCE
%
% Derives snow-thickness-dependent sea-ice density parameterizations shown in Figure 2.
%
% The script compares two predictor choices:
%   1. Ice thickness and snow thickness.
%   2. Hydrostatic freeboard and snow thickness.
%
% For each case, the model is fitted using measured snow thickness, then only the snow-transition thresholds are refitted using ERA5 SnowModel-LG snow thickness.
%
% Input:
%   data/final/snow_model_matchups_with_models.mat
%
% Output:
%   figures/Fig2.png
%

close all

finalDir = fullfile(repoRoot, ...
    'data', 'final');

figureDir = fullfile(repoRoot, ...
    'figures');

colormapDir = fullfile(repoRoot, ...
    'data', 'colormaps');

if ~exist(figureDir, 'dir')
    mkdir(figureDir)
end

inputFile = fullfile(finalDir, ...
    'snow_model_matchups_with_models.mat');

load(inputFile,'D')

S = load(fullfile(colormapDir,'lipari.mat'));
lipari = S.lipari;

densityMode = 1;
% Density definition:
% 1 = laboratory measurements
% 2 = in-situ connected porosity
% 3 = in-situ disconnected porosity

rho_cols = {"density_lab_kgm3", ...
            "density_insitu_connected_kgm3", ...
            "density_insitu_disconnected_kgm3"};
rho_names = {"lab","connected","disconnected"};
rho_col = rho_cols{densityMode};

hs_color_max = 0.08;
hi_lim = [0 3.5];
hf_lim = [0 0.35];
rho_lim = [850 930];
rho_w = 1024; % Seawater density used for hydrostatic freeboard estimation (kg m^-3).

% Enhance color contrast at low snow thickness.
ntrim = 25;
lipari = lipari(1:end-ntrim,:);
lipari = flipud(lipari);
lipari = interp1(linspace(0,1,size(lipari,1)),lipari,linspace(0,1,256));

requiredVars = ["date","ice_age","dataset","ice_thickness_m","hs","hs_smlg",rho_col];
missingVars = requiredVars(~ismember(requiredVars,string(D.Properties.VariableNames)));
if ~isempty(missingVars)
    error('D is missing required variables: %s',strjoin(cellstr(missingVars),', '))
end

dataset_all = string(D.dataset);

% Exclude Svalbard fjord observations from the snow-dependence analysis because this parameterization targets Arctic pack-ice conditions.
is_svalbard = contains(lower(dataset_all),"svalbard");
mask = ~is_svalbard;

hi_all = D.ice_thickness_m(mask);
rho_all = D.(rho_col)(mask);
date_all = D.date(mask);
age_all = D.ice_age(mask);
hs_meas0 = D.hs(mask);
hs_smlg0 = D.hs_smlg(mask);

% Density parameterization as a function of ice thickness and snow thickness.
[x_hi_meas,hs_hi_meas,rho_hi_meas,valid_hi_meas] = hi_prepare_thickness_data(hi_all,rho_all,hs_meas0);
[x_hi_smlg,hs_hi_smlg,rho_hi_smlg,valid_hi_smlg] = hi_prepare_thickness_data(hi_all,rho_all,hs_smlg0);

age_hi_meas = age_all(valid_hi_meas);
age_hi_smlg = age_all(valid_hi_smlg);

hi_ref_thick = 1.1;
rho_thick_at_ref = 911;
rho_thick_high_manual = 908.5;
hi_thick_high = 3.5;
rho_thick_0_bounds = [890 915];
hi0_low_bounds = [1.2 2.3];
rho_low_const_bounds_hi = [890 910];

% Bounds for the snow-thickness transition zone between low-snow and high-snow density branches.
hs_low_min = 0.00;
hs_low_max = 0.03;
hs_high_min = 0.05;
hs_high_max = 0.12;

p0_hi = [857, 908, 0.015, 0.07, 1.9, 905];
% Bootstrap resampling parameters.
nboot_hi = 1000;
rng(1)

% Fit the ice-thickness-based density model using measured snow thickness.
[p_hi_meas,RMSE_hi_meas,R2_hi_meas] = hi_fit_panel1_model( ...
    x_hi_meas,hs_hi_meas,rho_hi_meas,p0_hi, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
    rho_thick_0_bounds,hi0_low_bounds,rho_low_const_bounds_hi, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

% Estimate parameter uncertainty using bootstrap resampling.
B_hi_meas = hi_bootstrap_panel1( ...
    x_hi_meas,hs_hi_meas,rho_hi_meas,nboot_hi,p0_hi, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
    rho_thick_0_bounds,hi0_low_bounds,rho_low_const_bounds_hi, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

% Refit snow-thickness transition thresholds using ERA5 SM-LG snow thickness.
[p_hi_smlg,RMSE_hi_smlg,R2_hi_smlg] = hi_fit_panel2_thresholds_only( ...
    x_hi_smlg,hs_hi_smlg,rho_hi_smlg,p_hi_meas, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

% Estimate uncertainty in the ERA5 SM-LG threshold parameters.
B_hi_smlg = hi_bootstrap_panel2_thresholds_only( ...
    x_hi_smlg,hs_hi_smlg,rho_hi_smlg,nboot_hi,p_hi_meas,p_hi_smlg, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

stats_hi_meas = hi_bootstrap_skill_ci(B_hi_meas,x_hi_meas,hs_hi_meas,rho_hi_meas,RMSE_hi_meas,R2_hi_meas, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual);
stats_hi_smlg = hi_bootstrap_skill_ci_with_ref(B_hi_smlg,B_hi_meas,x_hi_smlg,hs_hi_smlg,rho_hi_smlg,RMSE_hi_smlg,R2_hi_smlg, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual);

xline_hi_meas = linspace(hi_lim(1),hi_lim(2),300);
hs_cases_hi_meas = [p_hi_meas(3), mean(p_hi_meas(3:4)), p_hi_meas(4)];
xline_hi_smlg = linspace(hi_lim(1),hi_lim(2),300);
hs_cases_hi_smlg = [p_hi_smlg(3), mean(p_hi_smlg(3:4)), p_hi_smlg(4)];

% Density parameterization as a function of freeboard and snow thickness.
[x_hf_meas,hs_hf_meas,rho_hf_meas,valid_hf_meas] = hf_prepare_freeboard_data(hi_all,rho_all,hs_meas0,date_all,rho_w);
[x_hf_smlg,hs_hf_smlg,rho_hf_smlg,valid_hf_smlg] = hf_prepare_freeboard_data(hi_all,rho_all,hs_smlg0,date_all,rho_w);

age_hf_meas = age_all(valid_hf_meas);
age_hf_smlg = age_all(valid_hf_smlg);

hf0_low_bounds = [0.15 0.35];
rho_low_const_bounds_hf = [895 905];
hf_ref_thick = 0.20;
rho_thick_const = 911;
hf_thick_high = 0.35;
rho_thick_high_manual_hf = 908;

p0_hf = [857, 908, 0.015, 0.07, 0.23];
% Bootstrap resampling parameters.
nboot_hf = 1000;
rng(1)

% Fit the freeboard-based density model using measured snow thickness.
[p_hf_meas,RMSE_hf_meas,R2_hf_meas] = hf_fit_full_model( ...
    x_hf_meas,hs_hf_meas,rho_hf_meas,p0_hf, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf, ...
    hf0_low_bounds,rho_low_const_bounds_hf, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

% Estimate parameter uncertainty using bootstrap resampling.
B_hf_meas = hf_bootstrap_full_model( ...
    x_hf_meas,hs_hf_meas,rho_hf_meas,nboot_hf,p0_hf, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf, ...
    hf0_low_bounds,rho_low_const_bounds_hf, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

% Refit snow-thickness transition thresholds using ERA5 SM-LG snow thickness.
[p_hf_smlg,RMSE_hf_smlg,R2_hf_smlg] = hf_fit_thresholds_only( ...
    x_hf_smlg,hs_hf_smlg,rho_hf_smlg,p_hf_meas, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

% Estimate uncertainty in the ERA5 SM-LG threshold parameters.
B_hf_smlg = hf_bootstrap_thresholds_only( ...
    x_hf_smlg,hs_hf_smlg,rho_hf_smlg,nboot_hf,p_hf_meas,p_hf_smlg, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

stats_hf_meas = hf_get_stats(B_hf_meas,RMSE_hf_meas,R2_hf_meas,numel(rho_hf_meas));
stats_hf_smlg = hf_get_stats(B_hf_smlg,RMSE_hf_smlg,R2_hf_smlg,numel(rho_hf_smlg));

xline_hf_meas = linspace(hf_lim(1),hf_lim(2),300);
hs_cases_hf_meas = [p_hf_meas(3), mean(p_hf_meas(3:4)), p_hf_meas(4)];
xline_hf_smlg = linspace(hf_lim(1),hf_lim(2),300);
hs_cases_hf_smlg = [p_hf_smlg(3), mean(p_hf_smlg(3:4)), p_hf_smlg(4)];

% Generate the four-panel comparison figure.
fig = figure;
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
hi_plot_panel(x_hi_meas,rho_hi_meas,hs_hi_meas,age_hi_meas,p_hi_meas,B_hi_meas,[],xline_hi_meas,hs_cases_hi_meas, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
    lipari,hs_color_max, ...
    sprintf('Measured snow thickness: R^2 = %.2f, RMSE = %.1f kg m^{-3}',R2_hi_meas,RMSE_hi_meas), ...
    stats_hi_meas,hi_lim,rho_lim,"full");
xlabel('Sea-ice thickness{\it h_i} (m)','Interpreter','tex','FontSize',11)

ax2 = nexttile;
hi_plot_panel(x_hi_smlg,rho_hi_smlg,hs_hi_smlg,age_hi_smlg,p_hi_smlg,B_hi_smlg,B_hi_meas,xline_hi_smlg,hs_cases_hi_smlg, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
    lipari,hs_color_max, ...
    sprintf('ERA5 SM-LG snow thickness: R^2 = %.2f, RMSE = %.1f kg m^{-3}',R2_hi_smlg,RMSE_hi_smlg), ...
    stats_hi_smlg,hi_lim,rho_lim,"thresholds");
xlabel('Sea-ice thickness{\it h_i} (m)','Interpreter','tex','FontSize',11)

ax3 = nexttile;
hf_plot_panel(x_hf_meas,rho_hf_meas,hs_hf_meas,age_hf_meas,p_hf_meas,B_hf_meas,[],xline_hf_meas,hs_cases_hf_meas, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf, ...
    lipari,hs_color_max, ...
    sprintf('Measured snow thickness: R^2 = %.2f, RMSE = %.1f kg m^{-3}',R2_hf_meas,RMSE_hf_meas), ...
    stats_hf_meas,hf_lim,rho_lim,"full");
xlabel('Sea-ice freeboard{\it h_f} (m)','Interpreter','tex','FontSize',11)

ax4 = nexttile;
hf_plot_panel(x_hf_smlg,rho_hf_smlg,hs_hf_smlg,age_hf_smlg,p_hf_smlg,B_hf_smlg,B_hf_meas,xline_hf_smlg,hs_cases_hf_smlg, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf, ...
    lipari,hs_color_max, ...
    sprintf('ERA5 SM-LG snow thickness: R^2 = %.2f, RMSE = %.1f kg m^{-3}',R2_hf_smlg,RMSE_hf_smlg), ...
    stats_hf_smlg,hf_lim,rho_lim,"thresholds");
xlabel('Sea-ice freeboard{\it h_f} (m)','Interpreter','tex','FontSize',11)

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Snow thickness {\it h_s} (m)';
cb.Label.Interpreter = 'tex';
cb.Limits = [0 hs_color_max];

fsz = 11;
text(ax1,-0.125,1.035,'(a)','Units','normalized','FontSize',fsz,'FontWeight','normal')
text(ax2,-0.125,1.035,'(b)','Units','normalized','FontSize',fsz,'FontWeight','normal')
text(ax3,-0.125,1.035,'(c)','Units','normalized','FontSize',fsz,'FontWeight','normal')
text(ax4,-0.125,1.035,'(d)','Units','normalized','FontSize',fsz,'FontWeight','normal')

set(fig,'Units','inches','Position',[1 1 12.0 10.0])

outputFigure = fullfile(figureDir, ...
    'Fig2.png');

set(findall(fig,'Type','axes'),'Toolbar',[])

exportgraphics(fig, outputFigure, 'Resolution',300)

close(fig)

fprintf('Generated Figure 2 using %s density.\n', ...
    rho_names{densityMode})

fprintf('Saved figure to:\n%s\n', ...
    outputFigure)

%% Helpers

function [x,hs,rho,valid] = hi_prepare_thickness_data(hi,rho,hs)
% Remove invalid observations and format ice-thickness model inputs.

x = hi(:);
hs = hs(:);
rho = rho(:);

valid = isfinite(x) & isfinite(hs) & isfinite(rho);

x = x(valid);
hs = hs(valid);
rho = rho(valid);

end

function [p_opt,RMSE,R2] = hi_fit_panel1_model( ...
    x,hs,rho,p0, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
    rho_thick_0_bounds,hi0_low_bounds,rho_low_const_bounds, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max)
% Fit the full ice-thickness/snow-thickness density model.

opts = optimset('MaxFunEvals',5e4,'MaxIter',5e4,'Display','off');

cost_fun = @(p) sqrt(mean((rho - hi_rho_model( ...
    p,x,hs,hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual)).^2)) + ...
    hi_penalty_panel1(p,hi0_low_bounds,rho_low_const_bounds,rho_thick_0_bounds, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

p_opt = fminsearch(cost_fun,p0,opts);

p_opt(2) = min(max(p_opt(2),rho_low_const_bounds(1)),rho_low_const_bounds(2));
p_opt(3) = min(max(p_opt(3),hs_low_min),hs_low_max);
p_opt(4) = min(max(p_opt(4),hs_high_min),hs_high_max);
p_opt(5) = min(max(p_opt(5),hi0_low_bounds(1)),hi0_low_bounds(2));
p_opt(6) = min(max(p_opt(6),rho_thick_0_bounds(1)),rho_thick_0_bounds(2));

rho_fit = hi_rho_model( ...
    p_opt,x,hs,hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual);

RMSE = sqrt(mean((rho - rho_fit).^2));
R2 = 1 - sum((rho - rho_fit).^2) ./ sum((rho - mean(rho)).^2);

end

function [p_opt,RMSE,R2] = hi_fit_panel2_thresholds_only( ...
    x,hs,rho,p_fixed, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max)
% Refit only snow-thickness transition thresholds for SM-LG snow input.

opts = optimset('MaxFunEvals',5e4,'MaxIter',5e4,'Display','off');

q0 = p_fixed(3:4);

cost_fun = @(q) sqrt(mean((rho - hi_rho_model( ...
    hi_replace_thresholds(p_fixed,q),x,hs,hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual)).^2)) + ...
    hi_penalty_thresholds(q,hs_low_min,hs_low_max,hs_high_min,hs_high_max);

q_opt = fminsearch(cost_fun,q0,opts);

q_opt(1) = min(max(q_opt(1),hs_low_min),hs_low_max);
q_opt(2) = min(max(q_opt(2),hs_high_min),hs_high_max);

if q_opt(1) > q_opt(2) - 0.005
    q_opt(2) = min(max(q_opt(1) + 0.005,hs_high_min),hs_high_max);
    q_opt(1) = min(q_opt(1),q_opt(2) - 0.005);
    q_opt(1) = min(max(q_opt(1),hs_low_min),hs_low_max);
end

p_opt = hi_replace_thresholds(p_fixed,q_opt);

rho_fit = hi_rho_model( ...
    p_opt,x,hs,hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual);

RMSE = sqrt(mean((rho - rho_fit).^2));
R2 = 1 - sum((rho - rho_fit).^2) ./ sum((rho - mean(rho)).^2);

end

function p = hi_replace_thresholds(p,q)

p = p(:).';
p(3) = q(1);
p(4) = q(2);

end

function rho_m = hi_rho_model( ...
    p,x,hs,hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual)
% Evaluate the ice-thickness/snow-thickness density parameterization.

rho_low_0 = p(1);
rho_low_const = p(2);
hs_low = p(3);
hs_high = p(4);
hi0_low = p(5);
rho_thick_0 = p(6);

rho_low = rho_low_0 + ...
    (rho_low_const - rho_low_0) .* min(max(x ./ hi0_low,0),1);

rho_thick = nan(size(x));

idx1 = x <= hi_ref_thick;
rho_thick(idx1) = rho_thick_0 + ...
    (rho_thick_at_ref - rho_thick_0) .* x(idx1) ./ hi_ref_thick;

idx2 = x > hi_ref_thick;
rho_thick(idx2) = rho_thick_at_ref + ...
    (rho_thick_high_manual - rho_thick_at_ref) .* ...
    (x(idx2) - hi_ref_thick) ./ (hi_thick_high - hi_ref_thick);

denom = max(hs_high - hs_low,1e-6);
w = min(max((hs - hs_low) ./ denom,0),1);

rho_m = (1 - w).*rho_low + w.*rho_thick;

end

function penalty = hi_penalty_panel1(p,hi0_low_bounds,rho_low_const_bounds,rho_thick_0_bounds, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max)
% Penalize parameter values outside physically plausible bounds during unconstrained fminsearch optimization.

rho_low_0 = p(1);
rho_low_const = p(2);
hs_low = p(3);
hs_high = p(4);
hi0_low = p(5);
rho_thick_0 = p(6);

penalty = 0;

penalty = penalty + 1e6*max(0,830-rho_low_0).^2;
penalty = penalty + 1e6*max(0,rho_low_0-900).^2;

penalty = penalty + 1e6*max(0,rho_low_const_bounds(1)-rho_low_const).^2;
penalty = penalty + 1e6*max(0,rho_low_const-rho_low_const_bounds(2)).^2;

penalty = penalty + hi_penalty_thresholds([hs_low hs_high], ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

penalty = penalty + 1e6*max(0,hi0_low_bounds(1)-hi0_low).^2;
penalty = penalty + 1e6*max(0,hi0_low-hi0_low_bounds(2)).^2;

penalty = penalty + 1e6*max(0,rho_thick_0_bounds(1)-rho_thick_0).^2;
penalty = penalty + 1e6*max(0,rho_thick_0-rho_thick_0_bounds(2)).^2;

end

function penalty = hi_penalty_thresholds(q,hs_low_min,hs_low_max,hs_high_min,hs_high_max)
% Penalize parameter values outside physically plausible bounds during unconstrained fminsearch optimization.

hs_low = q(1);
hs_high = q(2);

penalty = 0;
penalty = penalty + 1e6*max(0,hs_low_min-hs_low).^2;
penalty = penalty + 1e6*max(0,hs_low-hs_low_max).^2;
penalty = penalty + 1e6*max(0,hs_high_min-hs_high).^2;
penalty = penalty + 1e6*max(0,hs_high-hs_high_max).^2;
penalty = penalty + 1e6*max(0,hs_low-hs_high+0.005).^2;

end

function B = hi_bootstrap_panel1( ...
    x,hs,rho,nboot,p0, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
    rho_thick_0_bounds,hi0_low_bounds,rho_low_const_bounds, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max)
% Bootstrap uncertainty for the full ice-thickness/snow-thickness model.

n = numel(rho);
B = nan(nboot,6);

for b = 1:nboot

    idx = randi(n,n,1);

    [p_b,~,~] = hi_fit_panel1_model( ...
        x(idx),hs(idx),rho(idx),p0, ...
        hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
        rho_thick_0_bounds,hi0_low_bounds,rho_low_const_bounds, ...
        hs_low_min,hs_low_max,hs_high_min,hs_high_max);

    if any(~isfinite(p_b)) || p_b(3) >= p_b(4)
        continue
    end

    B(b,:) = p_b;

end

B = B(all(isfinite(B),2),:);

end

function B = hi_bootstrap_panel2_thresholds_only( ...
    x,hs,rho,nboot,p_fixed,p_start, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max)
% Bootstrap uncertainty when only SM-LG snow-threshold parameters are refit.

n = numel(rho);
B = nan(nboot,6);

p_fixed = p_fixed(:).';
p_start = p_start(:).';

for b = 1:nboot

    idx = randi(n,n,1);

    p_fixed_b = p_fixed;
    p_fixed_b(3:4) = p_start(3:4);

    [p_b,~,~] = hi_fit_panel2_thresholds_only( ...
        x(idx),hs(idx),rho(idx),p_fixed_b, ...
        hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
        hs_low_min,hs_low_max,hs_high_min,hs_high_max);

    if any(~isfinite(p_b)) || p_b(3) >= p_b(4)
        continue
    end

    B(b,:) = p_b;

end

B = B(all(isfinite(B),2),:);

end

function stats = hi_bootstrap_skill_ci(B,x,hs,rho,RMSE,R2, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual)

stats.RMSE = RMSE;
stats.R2 = R2;
stats.N = numel(rho);

if isempty(B)
    stats.RMSE_ci = [nan nan];
    stats.R2_ci = [nan nan];
    stats.hs_low_ci = [nan nan];
    stats.hs_high_ci = [nan nan];
    return
end

RMSE_b = nan(size(B,1),1);
R2_b = nan(size(B,1),1);

for b = 1:size(B,1)

    rho_fit = hi_rho_model( ...
        B(b,:),x,hs,hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual);

    RMSE_b(b) = sqrt(mean((rho - rho_fit).^2));
    R2_b(b) = 1 - sum((rho - rho_fit).^2) ./ sum((rho - mean(rho)).^2);

end

stats.RMSE_ci = prctile(RMSE_b,[2.5 97.5]);
stats.R2_ci = prctile(R2_b,[2.5 97.5]);
stats.hs_low_ci = prctile(B(:,3),[2.5 97.5]);
stats.hs_high_ci = prctile(B(:,4),[2.5 97.5]);

end

function stats = hi_bootstrap_skill_ci_with_ref(B,B_ref,x,hs,rho,RMSE,R2, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual)

stats.RMSE = RMSE;
stats.R2 = R2;
stats.N = numel(rho);

if isempty(B)
    stats.RMSE_ci = [nan nan];
    stats.R2_ci = [nan nan];
    stats.hs_low_ci = [nan nan];
    stats.hs_high_ci = [nan nan];
    return
end

RMSE_b = nan(size(B,1),1);
R2_b = nan(size(B,1),1);

for b = 1:size(B,1)

    p_b = B(b,:);

    if ~isempty(B_ref)
        jj = randi(size(B_ref,1));
        p_b([1 2 5 6]) = B_ref(jj,[1 2 5 6]);
    end

    rho_fit = hi_rho_model( ...
        p_b,x,hs,hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual);

    RMSE_b(b) = sqrt(mean((rho - rho_fit).^2));
    R2_b(b) = 1 - sum((rho - rho_fit).^2) ./ sum((rho - mean(rho)).^2);

end

stats.RMSE_ci = prctile(RMSE_b,[2.5 97.5]);
stats.R2_ci = prctile(R2_b,[2.5 97.5]);
stats.hs_low_ci = prctile(B(:,3),[2.5 97.5]);
stats.hs_high_ci = prctile(B(:,4),[2.5 97.5]);

end

function hi_plot_panel(x,rho,hs,age,p,B,B_ref,x_line,hs_cases, ...
    hi_ref_thick,rho_thick_at_ref,hi_thick_high,rho_thick_high_manual, ...
    lipari,hs_color_max,title_text,stats,x_lim,rho_lim,text_mode)

hold on

isFYI = strcmpi(string(age),'FYI');
isSYI = strcmpi(string(age),'SYI') | strcmpi(string(age),'MYI');

scatter(x(isFYI),rho(isFYI),40,min(hs(isFYI),hs_color_max), ...
    'filled','MarkerEdgeColor',[.5 .5 .5],'LineWidth',0.7,'HandleVisibility','off');

scatter(x(isSYI),rho(isSYI),40,min(hs(isSYI),hs_color_max), ...
    'o','MarkerFaceColor','none','MarkerEdgeColor','flat','LineWidth',1.0,'HandleVisibility','off');

hFYI = scatter(nan,nan,40,'k','filled','DisplayName','FYI');
hSYI = scatter(nan,nan,40,'k','o','MarkerFaceColor','none','LineWidth',1.0,'DisplayName','SYI & MYI');

colormap(gca,lipari)
clim([0 hs_color_max])

hs_labels = {sprintf('h_s = %.2f m',round(max(hs_cases(1),0),2)), ...
             sprintf('h_s = %.2f m',round(max(hs_cases(2),0),2)), ...
             sprintf('h_s = %.2f m',round(max(hs_cases(3),0),2))};

h_lines = gobjects(numel(hs_cases),1);
h_band = gobjects(1,1);

for i = 1:numel(hs_cases)

    line_color = hi_get_colormap_color(lipari,hs_cases(i),0,hs_color_max);
    rho_boot = nan(numel(x_line),size(B,1));

    for b = 1:size(B,1)

        p_b = B(b,:);

        if ~isempty(B_ref)
            jj = randi(size(B_ref,1));
            p_b([1 2 5 6]) = B_ref(jj,[1 2 5 6]);
        end

        if i == 1
            hs_b = p_b(3);
        elseif i == 2
            hs_b = mean(p_b(3:4));
        else
            hs_b = p_b(4);
        end

        rho_boot(:,b) = hi_rho_model( ...
            p_b,x_line,hs_b,hi_ref_thick,rho_thick_at_ref, ...
            hi_thick_high,rho_thick_high_manual);
    end

    rho_lo = prctile(rho_boot,2.5,2);
    rho_hi = prctile(rho_boot,97.5,2);

    hb = fill([x_line fliplr(x_line)], ...
        [rho_lo' fliplr(rho_hi')], ...
        [0.75 0.75 0.75], ...
        'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');

    if i == 1
        h_band = hb;
    end

    rho_line = hi_rho_model( ...
        p,x_line,hs_cases(i),hi_ref_thick,rho_thick_at_ref, ...
        hi_thick_high,rho_thick_high_manual);

    h_lines(i) = plot(x_line,rho_line,'LineWidth',2.5, ...
        'Color',line_color,'DisplayName',hs_labels{i});
end

ylabel('Sea-ice density \rho (kg m^{-3})','Interpreter','tex','FontSize',11)
title(title_text,'Interpreter','tex','FontWeight','normal','FontSize',11)

if text_mode == "full"
    hi_add_model_text_full(gca,p,hi_ref_thick,rho_thick_at_ref, ...
        hi_thick_high,rho_thick_high_manual,stats)
elseif text_mode == "thresholds"
    hi_add_model_text_thresholds(gca,p,stats)
end

legend([h_lines; h_band; hFYI; hSYI], ...
    [hs_labels(:); {'95% CI'; 'FYI'; 'SYI & MYI'}], ...
    'Location','north','NumColumns',3,'Box','off')

xlim(x_lim)
ylim(rho_lim)
box on

end

function hi_add_model_text_full(ax,p,hi_ref_thick,rho_thick_at_ref, ...
    hi_thick_high,rho_thick_high_manual,stats)

rho_low_0 = p(1);
rho_low_const = p(2);
hs_low = p(3);
hs_high = p(4);
hi0_low = p(5);
rho_thick_0 = p(6);

b_low = (rho_low_const - rho_low_0) ./ hi0_low;
b_thick = (rho_thick_high_manual - rho_thick_at_ref) ./ (hi_thick_high - hi_ref_thick);

hs_tol = 0.005;

if abs(hs_low) < hs_tol
    wtxt = sprintf([ ...
        '$w=\\frac{h_s}{%.2f},\\quad 0<h_s<%.2f\\,\\mathrm{m}$\n' ...
        '$w=1,\\quad h_s\\geq %.2f\\,\\mathrm{m}$\n'], ...
        hs_high,hs_high,hs_high);
else
    wtxt = sprintf([ ...
        '$w=0,\\quad h_s\\leq %.2f\\,\\mathrm{m}$\n' ...
        '$w=\\frac{h_s-%.2f}{%.2f-%.2f},\\quad %.2f<h_s<%.2f\\,\\mathrm{m}$\n' ...
        '$w=1,\\quad h_s\\geq %.2f\\,\\mathrm{m}$\n'], ...
        hs_low, ...
        hs_low,hs_high,hs_low,hs_low,hs_high, ...
        hs_high);
end

txt = sprintf([ ...
    '$\\rho(h_i,h_s)=(1-w)\\rho_{\\ell}(h_i)+w\\rho_h(h_i)$\n' ...
    '$\\rho_{\\ell}=%.0f%+.0f h_i,\\quad h_i<%.1f\\,\\mathrm{m}$\n' ...
    '$\\rho_{\\ell}=%.0f,\\quad h_i\\geq %.1f\\,\\mathrm{m}$\n' ...
    '$\\rho_h=%.0f+\\frac{%.0f-%.0f}{%.1f}h_i,\\quad h_i<%.1f\\,\\mathrm{m}$\n' ...
    '$\\rho_h=%.0f%+.1f\\,(h_i-%.1f),\\quad h_i\\geq%.1f\\,\\mathrm{m}$\n' ...
    '%s' ...
    '$R^2=%.2f,\\; RMSE=%.1f\\,\\mathrm{kg\\,m^{-3}},\\; N=%d$'], ...
    rho_low_0,b_low,hi0_low, ...
    rho_low_const,hi0_low, ...
    rho_thick_0,rho_thick_at_ref,rho_thick_0,hi_ref_thick,hi_ref_thick, ...
    rho_thick_at_ref,b_thick,hi_ref_thick,hi_ref_thick, ...
    wtxt, ...
    stats.R2, ...
    stats.RMSE,stats.N);

text(ax,0.98,0.01,txt, ...
    'Units','normalized','Interpreter','latex', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'FontSize',10);

end

function hi_add_model_text_thresholds(ax,p,stats)

hs_low = p(3);
hs_high = p(4);

wtxt = make_w_text(hs_low,hs_high);

txt = sprintf([ ...
    '$\\rho_{\\ell}(h_i),\\rho_h(h_i)\\;\\mathrm{fixed\\;from\\;panel\\;(a)}$\n' ...
    '%s' ...
    '$R^2=%.2f,\\; RMSE=%.1f\\,\\mathrm{kg\\,m^{-3}},\\; N=%d$'], ...
    wtxt, ...
    stats.R2, ...
    stats.RMSE,stats.N);

text(ax,0.98,0.01,txt, ...
    'Units','normalized','Interpreter','latex', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'FontSize',10);

end

function c = hi_get_colormap_color(cmap,value,cmin,cmax)

value = min(max(value,cmin),cmax);

if cmax == cmin
    idx = 1;
else
    idx = 1 + round((value - cmin) ./ (cmax - cmin) .* (size(cmap,1) - 1));
end

idx = min(max(idx,1),size(cmap,1));
c = cmap(idx,:);

end

function [x,hs,rho,valid] = hf_prepare_freeboard_data(hi,rho,hs,date,rho_w)
% Estimate hydrostatic freeboard and prepare freeboard-model inputs.

aug1_this_year = datetime(year(date),8,1);
aug1_next = aug1_this_year;
aug1_next(date > aug1_this_year) = datetime(year(date(date > aug1_this_year)) + 1,8,1);

% Seasonal snow-density parameterization (Mallett, 2025) used for hydrostatic freeboard estimation.
rho_s = 0.35 .* days(aug1_next - date) + 239.78;

x = hi - (rho.*hi + rho_s.*hs)./rho_w;
x = max(x,0);

x = x(:);
hs = hs(:);
rho = rho(:);

valid = isfinite(x) & isfinite(hs) & isfinite(rho);

x = x(valid);
hs = hs(valid);
rho = rho(valid);

end

function [p_opt,RMSE,R2] = hf_fit_full_model( ...
    x,hs,rho,p0, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual, ...
    hf0_low_bounds,rho_low_const_bounds, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max)

opts = optimset('MaxFunEvals',5e4,'MaxIter',5e4,'Display','off');

cost_fun = @(p) sqrt(mean((rho - hf_rho_model( ...
    p,x,hs,hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual)).^2)) + ...
    hf_penalty_full(p,hf0_low_bounds,rho_low_const_bounds, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

p_opt = fminsearch(cost_fun,p0,opts);

p_opt(2) = min(max(p_opt(2),rho_low_const_bounds(1)),rho_low_const_bounds(2));
p_opt(3) = min(max(p_opt(3),hs_low_min),hs_low_max);
p_opt(4) = min(max(p_opt(4),hs_high_min),hs_high_max);
p_opt(5) = min(max(p_opt(5),hf0_low_bounds(1)),hf0_low_bounds(2));

rho_fit = hf_rho_model(p_opt,x,hs,hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual);

RMSE = sqrt(mean((rho - rho_fit).^2));
R2 = 1 - sum((rho - rho_fit).^2) ./ sum((rho - mean(rho)).^2);

end

function [p_opt,RMSE,R2] = hf_fit_thresholds_only( ...
    x,hs,rho,p_fixed, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max)

opts = optimset('MaxFunEvals',5e4,'MaxIter',5e4,'Display','off');

q0 = p_fixed(3:4);

cost_fun = @(q) sqrt(mean((rho - hf_rho_model( ...
    hf_replace_thresholds(p_fixed,q),x,hs,hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual)).^2)) + ...
    hf_penalty_thresholds(q,hs_low_min,hs_low_max,hs_high_min,hs_high_max);

q_opt = fminsearch(cost_fun,q0,opts);
q_opt = hf_clip_thresholds(q_opt,hs_low_min,hs_low_max,hs_high_min,hs_high_max);

p_opt = hf_replace_thresholds(p_fixed,q_opt);

rho_fit = hf_rho_model(p_opt,x,hs,hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual);

RMSE = sqrt(mean((rho - rho_fit).^2));
R2 = 1 - sum((rho - rho_fit).^2) ./ sum((rho - mean(rho)).^2);

end

function p = hf_replace_thresholds(p,q)

p = p(:).';
p(3) = q(1);
p(4) = q(2);

end

function q = hf_clip_thresholds(q,hs_low_min,hs_low_max,hs_high_min,hs_high_max)

q(1) = min(max(q(1),hs_low_min),hs_low_max);
q(2) = min(max(q(2),hs_high_min),hs_high_max);

if q(1) > q(2) - 0.005
    q(2) = min(max(q(1) + 0.005,hs_high_min),hs_high_max);
    q(1) = min(q(1),q(2) - 0.005);
    q(1) = min(max(q(1),hs_low_min),hs_low_max);
end

end

function rho_m = hf_rho_model(p,x,hs,hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual)
% Evaluate the freeboard/snow-thickness density parameterization.

rho_low_0 = p(1);
rho_low_const = p(2);
hs_low = p(3);
hs_high = p(4);
hf0_low = p(5);

rho_low = rho_low_0 + ...
    (rho_low_const - rho_low_0) .* min(max(x ./ hf0_low,0),1);

rho_thick = rho_thick_const .* ones(size(x));

idx = x > hf_ref_thick;
rho_thick(idx) = rho_thick_const + ...
    (rho_thick_high_manual - rho_thick_const) .* ...
    min(max((x(idx) - hf_ref_thick) ./ (hf_thick_high - hf_ref_thick),0),1);

denom = max(hs_high - hs_low,1e-6);
w = min(max((hs - hs_low) ./ denom,0),1);

rho_m = (1 - w).*rho_low + w.*rho_thick;

end

function penalty = hf_penalty_full(p,hf0_low_bounds,rho_low_const_bounds, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max)

rho_low_0 = p(1);
rho_low_const = p(2);
hs_low = p(3);
hs_high = p(4);
hf0_low = p(5);

penalty = 0;

penalty = penalty + 1e6*max(0,830-rho_low_0).^2;
penalty = penalty + 1e6*max(0,rho_low_0-900).^2;

penalty = penalty + 1e6*max(0,rho_low_const_bounds(1)-rho_low_const).^2;
penalty = penalty + 1e6*max(0,rho_low_const-rho_low_const_bounds(2)).^2;

penalty = penalty + hf_penalty_thresholds([hs_low hs_high], ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

penalty = penalty + 1e6*max(0,hf0_low_bounds(1)-hf0_low).^2;
penalty = penalty + 1e6*max(0,hf0_low-hf0_low_bounds(2)).^2;

end

function penalty = hf_penalty_thresholds(q,hs_low_min,hs_low_max,hs_high_min,hs_high_max)

hs_low = q(1);
hs_high = q(2);

penalty = 0;
penalty = penalty + 1e6*max(0,hs_low_min-hs_low).^2;
penalty = penalty + 1e6*max(0,hs_low-hs_low_max).^2;
penalty = penalty + 1e6*max(0,hs_high_min-hs_high).^2;
penalty = penalty + 1e6*max(0,hs_high-hs_high_max).^2;
penalty = penalty + 1e6*max(0,hs_low-hs_high+0.005).^2;

end

function B = hf_bootstrap_full_model( ...
    x,hs,rho,nboot,p0, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual, ...
    hf0_low_bounds,rho_low_const_bounds, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max)
% Bootstrap uncertainty for the full freeboard/snow-thickness model.

n = numel(rho);
B = nan(nboot,5);

for b = 1:nboot

    idx = randi(n,n,1);

    [p_b,~,~] = hf_fit_full_model( ...
        x(idx),hs(idx),rho(idx),p0, ...
        hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual, ...
        hf0_low_bounds,rho_low_const_bounds, ...
        hs_low_min,hs_low_max,hs_high_min,hs_high_max);

    if any(~isfinite(p_b)) || p_b(3) >= p_b(4)
        continue
    end

    B(b,:) = p_b;

end

B = B(all(isfinite(B),2),:);

end

function B = hf_bootstrap_thresholds_only( ...
    x,hs,rho,nboot,p_fixed,p_start, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max)
% Bootstrap uncertainty when only SM-LG snow-threshold parameters are refit.

n = numel(rho);
B = nan(nboot,5);

p_fixed = p_fixed(:).';
p_start = p_start(:).';

for b = 1:nboot

    idx = randi(n,n,1);

    p_fixed_b = p_fixed;
    p_fixed_b(3:4) = p_start(3:4);

    [p_b,~,~] = hf_fit_thresholds_only( ...
        x(idx),hs(idx),rho(idx),p_fixed_b, ...
        hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual, ...
        hs_low_min,hs_low_max,hs_high_min,hs_high_max);

    if any(~isfinite(p_b)) || p_b(3) >= p_b(4)
        continue
    end

    B(b,:) = p_b;

end

B = B(all(isfinite(B),2),:);

end

function stats = hf_get_stats(B,RMSE,R2,N)

stats.RMSE = RMSE;
stats.R2 = R2;
stats.N = N;

if isempty(B)
    stats.hs_low_ci = [nan nan];
    stats.hs_high_ci = [nan nan];
else
    stats.hs_low_ci = prctile(B(:,3),[2.5 97.5]);
    stats.hs_high_ci = prctile(B(:,4),[2.5 97.5]);
end

end

function hf_plot_panel(x,rho,hs,age,p,B,B_ref,x_line,hs_cases, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual, ...
    lipari,hs_color_max,title_text,stats,x_lim,rho_lim,text_mode)

hold on

isFYI = strcmpi(string(age),'FYI');
isSYI = strcmpi(string(age),'SYI') | strcmpi(string(age),'MYI');

scatter(x(isFYI),rho(isFYI),40,min(hs(isFYI),hs_color_max), ...
    'filled','MarkerEdgeColor',[.5 .5 .5],'LineWidth',0.7,'HandleVisibility','off');

scatter(x(isSYI),rho(isSYI),40,min(hs(isSYI),hs_color_max), ...
    'o','MarkerFaceColor','none','MarkerEdgeColor','flat','LineWidth',1.0,'HandleVisibility','off');

hFYI = scatter(nan,nan,40,'k','filled','DisplayName','FYI');
hSYI = scatter(nan,nan,40,'k','o','MarkerFaceColor','none','LineWidth',1.0,'DisplayName','SYI & MYI');

colormap(gca,lipari)
clim([0 hs_color_max])

hs_labels = {sprintf('h_s = %.2f m',round(max(hs_cases(1),0),2)), ...
             sprintf('h_s = %.2f m',round(max(hs_cases(2),0),2)), ...
             sprintf('h_s = %.2f m',round(max(hs_cases(3),0),2))};

h_lines = gobjects(numel(hs_cases),1);
h_band = gobjects(1,1);

for i = 1:numel(hs_cases)

    line_color = hf_get_colormap_color(lipari,hs_cases(i),0,hs_color_max);
    rho_boot = nan(numel(x_line),size(B,1));

    for b = 1:size(B,1)

        p_b = B(b,:);

        if ~isempty(B_ref)
            jj = randi(size(B_ref,1));
            p_b([1 2 5]) = B_ref(jj,[1 2 5]);
        end

        if i == 1
            hs_b = p_b(3);
        elseif i == 2
            hs_b = mean(p_b(3:4));
        else
            hs_b = p_b(4);
        end

        rho_boot(:,b) = hf_rho_model( ...
            p_b,x_line,hs_b,hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual);
    end

    rho_lo = prctile(rho_boot,2.5,2);
    rho_hi = prctile(rho_boot,97.5,2);

    hb = fill([x_line fliplr(x_line)], ...
        [rho_lo' fliplr(rho_hi')], ...
        [0.75 0.75 0.75], ...
        'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');

    if i == 1
        h_band = hb;
    end

    rho_line = hf_rho_model( ...
        p,x_line,hs_cases(i),hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual);

    h_lines(i) = plot(x_line,rho_line,'LineWidth',2.5, ...
        'Color',line_color,'DisplayName',hs_labels{i});
end

ylabel('Sea-ice density \rho (kg m^{-3})','Interpreter','tex','FontSize',11)
title(title_text,'Interpreter','tex','FontWeight','normal','FontSize',11)

if text_mode == "full"
    hf_add_model_text_full(gca,p,hf_ref_thick,rho_thick_const, ...
        hf_thick_high,rho_thick_high_manual,stats)
elseif text_mode == "thresholds"
    hf_add_model_text_thresholds(gca,p,stats)
end

legend([h_lines; h_band; hFYI; hSYI], ...
    [hs_labels(:); {'95% CI'; 'FYI'; 'SYI & MYI'}], ...
    'Location','north','NumColumns',3,'Box','off')

xlim(x_lim)
ylim(rho_lim)
box on

end

function hf_add_model_text_full(ax,p,hf_ref_thick,rho_thick_const, ...
    hf_thick_high,rho_thick_high_manual,stats)

rho_low_0 = p(1);
rho_low_const = p(2);
hs_low = p(3);
hs_high = p(4);
hf0_low = p(5);

% Rounded display values for the low-snow branch.
rho_low_0_txt = round(rho_low_0);
rho_low_const_txt = round(rho_low_const);
hf0_low_txt = round(hf0_low,2);

b_low_txt = (rho_low_const_txt - rho_low_0_txt) ./ hf0_low_txt;

b_thick = (rho_thick_high_manual - rho_thick_const) ./ ...
    (hf_thick_high - hf_ref_thick);

wtxt = make_w_text(hs_low,hs_high);

txt = sprintf([ ...
    '$\\rho(h_f,h_s)=(1-w)\\rho_{\\ell}(h_f)+w\\rho_h(h_f)$\n' ...
    '$\\rho_{\\ell}=%.0f%+.0f h_f,\\quad h_f<%.2f\\,\\mathrm{m}$\n' ...
    '$\\rho_{\\ell}=%.0f,\\quad h_f\\geq %.2f\\,\\mathrm{m}$\n' ...
    '$\\rho_h=%.0f,\\quad h_f\\leq %.2f\\,\\mathrm{m}$\n' ...
    '$\\rho_h=%.0f%+.0f\\,(h_f-%.2f),\\quad h_f>%.2f\\,\\mathrm{m}$\n' ...
    '%s' ...
    '$R^2=%.2f,\\; RMSE=%.1f\\,\\mathrm{kg\\,m^{-3}},\\; N=%d$'], ...
    rho_low_0_txt,b_low_txt,hf0_low_txt, ...
    rho_low_const_txt,hf0_low_txt, ...
    rho_thick_const,hf_ref_thick, ...
    rho_thick_const,b_thick,hf_ref_thick,hf_ref_thick, ...
    wtxt, ...
    stats.R2,stats.RMSE,stats.N);

text(ax,0.98,0.01,txt, ...
    'Units','normalized','Interpreter','latex', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'FontSize',10);

end

function hf_add_model_text_thresholds(ax,p,stats)

hs_low = p(3);
hs_high = p(4);

wtxt = make_w_text(hs_low,hs_high);

txt = sprintf([ ...
    '$\\rho_{\\ell}(h_f),\\rho_h(h_f)\\;\\mathrm{fixed\\;from\\;panel\\;(c)}$\n' ...
    '%s' ...
    '$R^2=%.2f,\\; RMSE=%.1f\\,\\mathrm{kg\\,m^{-3}},\\; N=%d$'], ...
    wtxt, ...
    stats.R2,stats.RMSE,stats.N);

text(ax,0.98,0.01,txt, ...
    'Units','normalized','Interpreter','latex', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'FontSize',10);

end

function c = hf_get_colormap_color(cmap,value,cmin,cmax)

value = min(max(value,cmin),cmax);

if cmax == cmin
    idx = 1;
else
    idx = 1 + round((value - cmin) ./ (cmax - cmin) .* (size(cmap,1) - 1));
end

idx = min(max(idx,1),size(cmap,1));
c = cmap(idx,:);

end

function wtxt = make_w_text(hs_low,hs_high)

hs_tol = 0.005;

if abs(hs_low) < hs_tol
    wtxt = sprintf([ ...
        '$w=\\frac{h_s}{%.2f},\\quad 0<h_s<%.2f\\,\\mathrm{m}$\n' ...
        '$w=1,\\quad h_s\\geq %.2f\\,\\mathrm{m}$\n'], ...
        hs_high,hs_high,hs_high);
else
    wtxt = sprintf([ ...
        '$w=0,\\quad h_s\\leq %.2f\\,\\mathrm{m}$\n' ...
        '$w=\\frac{h_s-%.2f}{%.2f-%.2f},\\quad %.2f<h_s<%.2f\\,\\mathrm{m}$\n' ...
        '$w=1,\\quad h_s\\geq %.2f\\,\\mathrm{m}$\n'], ...
        hs_low, ...
        hs_low,hs_high,hs_low,hs_low,hs_high, ...
        hs_high);
end

end

end