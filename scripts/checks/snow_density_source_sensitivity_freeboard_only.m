%
% Snow-density source sensitivity for freeboard-based sea-ice density model.
%
% This standalone script tests exactly one question:
%
%   Does the performance of ERA5/MERRA2 snow-thickness parameterizations
%   change when their own SnowModel-LG snow density is used instead of the
%   simple seasonal snow-density parameterization?
%
% Controlled comparison:
%
%   ERA5:
%     (a) hs_smlg  + rhos_param
%     (b) hs_smlg  + rhos_smlg
%
%   MERRA2:
%     (c) hs_merra + rhos_param
%     (d) hs_merra + rhos_merra
%
% In each pair, snow thickness is kept fixed and only snow density changes.
%
% The script does not save or modify any .mat file. It only saves the figure.

clear
close all
clc

% USER SETTINGS

repoRoot = 'C:\Users\evsalg001\Documents\MATLAB\Density parametrization\script';

densityMode = 1; % 1 = lab, 2 = connected, 3 = disconnected

doBootstrap = false;   % true = slow but includes bootstrap CI bands
nboot_hf = 1000;

saveFigure = true;

outputFigureName = 'Fig2_snow_density_source_sensitivity.png';

rng(1)

% PATHS

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

outputFigure = fullfile(figureDir, ...
    outputFigureName);

% LOAD DATA

load(inputFile,'D')

S = load(fullfile(colormapDir,'lipari.mat'));
lipari = S.lipari;

rho_cols = {"density_lab_kgm3", ...
            "density_insitu_connected_kgm3", ...
            "density_insitu_disconnected_kgm3"};

rho_names = {"lab","connected","disconnected"};
rho_col = rho_cols{densityMode};

requiredVars = [ ...
    "date", ...
    "ice_age", ...
    "dataset", ...
    "ice_thickness_m", ...
    "hs", ...
    "hs_smlg", ...
    "hs_merra", ...
    "rhos_param", ...
    "rhos_smlg", ...
    "rhos_merra", ...
    rho_col];

missingVars = requiredVars(~ismember(requiredVars,string(D.Properties.VariableNames)));

if ~isempty(missingVars)
    error('D is missing required variables: %s',strjoin(cellstr(missingVars),', '))
end

if ~isdatetime(D.date)
    D.date = datetime(D.date);
end

% PLOT / MODEL SETTINGS

hs_color_max = 0.08;
hf_lim = [0 0.35];
rho_lim = [850 930];
rho_w = 1024; % kg m^-3

% Enhance color contrast at low snow thickness.
ntrim = 25;
lipari = lipari(1:end-ntrim,:);
lipari = flipud(lipari);
lipari = interp1(linspace(0,1,size(lipari,1)),lipari,linspace(0,1,256));

dataset_all = string(D.dataset);

if ismember("is_freeze",string(D.Properties.VariableNames))
    is_freeze = logical(D.is_freeze);
else
    is_freeze = false(height(D),1);
end

is_svalbard = contains(lower(dataset_all),"svalbard");

% Same filtering logic as the snow-dependence analysis.
mask = ~(is_freeze | is_svalbard);

hi_all = D.ice_thickness_m(mask);
rho_all = D.(rho_col)(mask);
age_all = D.ice_age(mask);

hs_meas0    = D.hs(mask);
hs_era50    = D.hs_smlg(mask);
hs_merra0   = D.hs_merra(mask);

rhos_param0 = D.rhos_param(mask);
rhos_era50  = D.rhos_smlg(mask);
rhos_merra0 = D.rhos_merra(mask);

% FREEBOARD MODEL SETTINGS

hf0_low_bounds = [0.15 0.35];
rho_low_const_bounds_hf = [895 905];

hf_ref_thick = 0.20;
rho_thick_const = 911;
hf_thick_high = 0.35;
rho_thick_high_manual_hf = 908;

hs_low_min = 0.00;
hs_low_max = 0.03;
hs_high_min = 0.05;
hs_high_max = 0.12;

p0_hf = [857, 908, 0.015, 0.07, 0.23];

% REFERENCE FIT: MEASURED SNOW THICKNESS + PARAMETERIZED SNOW DENSITY
%
% This calibrates the freeboard-density curves. The four comparison panels
% then keep these curves fixed and only refit the snow-thickness transition
% thresholds for each model-snow case.

[x_hf_ref,hs_hf_ref,rho_hf_ref,valid_hf_ref] = hf_prepare_freeboard_data( ...
    hi_all,rho_all,hs_meas0,rhos_param0,rho_w);

age_hf_ref = age_all(valid_hf_ref);

[p_hf_ref,RMSE_hf_ref,R2_hf_ref] = hf_fit_full_model( ...
    x_hf_ref,hs_hf_ref,rho_hf_ref,p0_hf, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf, ...
    hf0_low_bounds,rho_low_const_bounds_hf, ...
    hs_low_min,hs_low_max,hs_high_min,hs_high_max);

if doBootstrap
    B_hf_ref = hf_bootstrap_full_model( ...
        x_hf_ref,hs_hf_ref,rho_hf_ref,nboot_hf,p0_hf, ...
        hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf, ...
        hf0_low_bounds,rho_low_const_bounds_hf, ...
        hs_low_min,hs_low_max,hs_high_min,hs_high_max);
else
    B_hf_ref = p_hf_ref;
end

stats_hf_ref = hf_get_stats(B_hf_ref,RMSE_hf_ref,R2_hf_ref,numel(rho_hf_ref));

fprintf('\nReference freeboard fit: measured hs + parameterized rho_s\n')
hf_print_fit(p_hf_ref,stats_hf_ref,hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf)

% CONTROLLED MODEL CASES

caseNames = [ ...
    "ERA5 h_s + parameterized \rho_s", ...
    "ERA5 h_s + ERA5 \rho_s", ...
    "MERRA2 h_s + parameterized \rho_s", ...
    "MERRA2 h_s + MERRA2 \rho_s"];

caseShort = [ ...
    "ERA5_paramrho", ...
    "ERA5_modelrho", ...
    "MERRA2_paramrho", ...
    "MERRA2_modelrho"];

caseHs = { ...
    hs_era50, ...
    hs_era50, ...
    hs_merra0, ...
    hs_merra0};

caseRhos = { ...
    rhos_param0, ...
    rhos_era50, ...
    rhos_param0, ...
    rhos_merra0};

caseAge = cell(1,4);
caseX = cell(1,4);
caseHsClean = cell(1,4);
caseRho = cell(1,4);
caseP = cell(1,4);
caseB = cell(1,4);
caseStats = cell(1,4);
caseRMSE = nan(1,4);
caseR2 = nan(1,4);
caseN = nan(1,4);

for c = 1:4

    [x_c,hs_c,rho_c,valid_c] = hf_prepare_freeboard_data( ...
        hi_all,rho_all,caseHs{c},caseRhos{c},rho_w);

    age_c = age_all(valid_c);

    [p_c,RMSE_c,R2_c] = hf_fit_thresholds_only( ...
        x_c,hs_c,rho_c,p_hf_ref, ...
        hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf, ...
        hs_low_min,hs_low_max,hs_high_min,hs_high_max);

    if doBootstrap
        B_c = hf_bootstrap_thresholds_only( ...
            x_c,hs_c,rho_c,nboot_hf,p_hf_ref,p_c, ...
            hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf, ...
            hs_low_min,hs_low_max,hs_high_min,hs_high_max);
    else
        B_c = p_c;
    end

    stats_c = hf_get_stats(B_c,RMSE_c,R2_c,numel(rho_c));

    caseX{c} = x_c;
    caseHsClean{c} = hs_c;
    caseRho{c} = rho_c;
    caseAge{c} = age_c;
    caseP{c} = p_c;
    caseB{c} = B_c;
    caseStats{c} = stats_c;

    caseRMSE(c) = RMSE_c;
    caseR2(c) = R2_c;
    caseN(c) = numel(rho_c);

    fprintf('\n%s\n',caseNames(c))
    hf_print_fit(p_c,stats_c,hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf)
end

% PRINT DIRECT DENSITY-SOURCE EFFECT

fprintf('\n============================================================\n')
fprintf('Controlled effect of replacing parameterized rho_s with model rho_s\n')
fprintf('Same snow thickness is used within each pair.\n')
fprintf('============================================================\n')

fprintf('\nERA5:\n')
fprintf('  hs_smlg + rhos_param: R2 = %.3f, RMSE = %.2f kg m-3, N = %d\n', ...
    caseR2(1),caseRMSE(1),caseN(1))
fprintf('  hs_smlg + rhos_smlg:  R2 = %.3f, RMSE = %.2f kg m-3, N = %d\n', ...
    caseR2(2),caseRMSE(2),caseN(2))
fprintf('  Change with ERA5 rho_s: dR2 = %.3f, dRMSE = %.2f kg m-3\n', ...
    caseR2(2)-caseR2(1),caseRMSE(2)-caseRMSE(1))

fprintf('\nMERRA2:\n')
fprintf('  hs_merra + rhos_param: R2 = %.3f, RMSE = %.2f kg m-3, N = %d\n', ...
    caseR2(3),caseRMSE(3),caseN(3))
fprintf('  hs_merra + rhos_merra: R2 = %.3f, RMSE = %.2f kg m-3, N = %d\n', ...
    caseR2(4),caseRMSE(4),caseN(4))
fprintf('  Change with MERRA2 rho_s: dR2 = %.3f, dRMSE = %.2f kg m-3\n', ...
    caseR2(4)-caseR2(3),caseRMSE(4)-caseRMSE(3))

%% FIGURE

xline_hf = linspace(hf_lim(1),hf_lim(2),300);

fig = figure;
set(fig,'Units','inches','Position',[1 1 12.0 8.8]);

tl = tiledlayout(2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

ax = gobjects(4,1);

for c = 1:4

    ax(c) = nexttile;
    p_c = caseP{c};
    hs_cases = [p_c(3), mean(p_c(3:4)), p_c(4)];

    hf_plot_panel( ...
        caseX{c},caseRho{c},caseHsClean{c},caseAge{c}, ...
        p_c,caseB{c},B_hf_ref,xline_hf,hs_cases, ...
        hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual_hf, ...
        lipari,hs_color_max, ...
        sprintf('%s: R^2 = %.2f, RMSE = %.1f kg m^{-3}', ...
        caseNames(c),caseR2(c),caseRMSE(c)), ...
        caseStats{c},hf_lim,rho_lim,"thresholds",doBootstrap);

    xlabel('Sea-ice freeboard {\it h_f} (m)', ...
        'Interpreter','tex', ...
        'FontSize',11)

    text(ax(c),-0.13,1.035,sprintf('(%c)','a'+c-1), ...
        'Units','normalized', ...
        'FontSize',11, ...
        'FontWeight','normal')
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Snow thickness {\it h_s} (m)';
cb.Label.Interpreter = 'tex';
cb.Limits = [0 hs_color_max];

set(findall(fig,'Type','axes'),'Toolbar',[])

if saveFigure
    exportgraphics(fig, outputFigure, 'Resolution',300)
    fprintf('\nSaved figure to:\n%s\n', outputFigure)
end

fprintf('\nGenerated snow-density source sensitivity figure using %s density.\n', ...
    rho_names{densityMode})

%% ========================================================================
% Helpers
% ========================================================================

function [x,hs,rho,valid] = hf_prepare_freeboard_data(hi,rho,hs,rho_s,rho_w)
% Estimate hydrostatic freeboard using explicitly supplied snow density.
%
% This is the central helper for this sensitivity test:
% snow thickness can be held fixed while snow density is changed.

hi = hi(:);
rho = rho(:);
hs = hs(:);
rho_s = rho_s(:);

x = hi - (rho.*hi + rho_s.*hs)./rho_w;
x = max(x,0);

valid = isfinite(x) & isfinite(hs) & isfinite(rho) & isfinite(rho_s);

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
    min(max((x(idx) - hf_ref_thick) ./ ...
    (hf_thick_high - hf_ref_thick),0),1);

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

if isempty(B) || size(B,1) <= 1
    stats.hs_low_ci = [nan nan];
    stats.hs_high_ci = [nan nan];
else
    stats.hs_low_ci = prctile(B(:,3),[2.5 97.5]);
    stats.hs_high_ci = prctile(B(:,4),[2.5 97.5]);
end

end

function hf_plot_panel(x,rho,hs,age,p,B,B_ref,x_line,hs_cases, ...
    hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual, ...
    lipari,hs_color_max,title_text,stats,x_lim,rho_lim,text_mode,doBootstrap)

hold on

isFYI = strcmpi(string(age),'FYI');
isSYI = strcmpi(string(age),'SYI') | strcmpi(string(age),'MYI');

scatter(x(isFYI),rho(isFYI),40,min(hs(isFYI),hs_color_max), ...
    'filled', ...
    'MarkerEdgeColor',[.5 .5 .5], ...
    'LineWidth',0.7, ...
    'HandleVisibility','off');

scatter(x(isSYI),rho(isSYI),40,min(hs(isSYI),hs_color_max), ...
    'o', ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor','flat', ...
    'LineWidth',1.0, ...
    'HandleVisibility','off');

hFYI = scatter(nan,nan,40,'k','filled','DisplayName','FYI');

hSYI = scatter(nan,nan,40,'k','o', ...
    'MarkerFaceColor','none', ...
    'LineWidth',1.0, ...
    'DisplayName','SYI & MYI');

colormap(gca,lipari)
clim([0 hs_color_max])

hs_labels = { ...
    sprintf('h_s = %.2f m',round(max(hs_cases(1),0),2)), ...
    sprintf('h_s = %.2f m',round(max(hs_cases(2),0),2)), ...
    sprintf('h_s = %.2f m',round(max(hs_cases(3),0),2))};

h_lines = gobjects(numel(hs_cases),1);
h_band = gobjects(1,1);

for i = 1:numel(hs_cases)

    line_color = hf_get_colormap_color(lipari,hs_cases(i),0,hs_color_max);

    if doBootstrap && size(B,1) > 1

        rho_boot = nan(numel(x_line),size(B,1));

        for b = 1:size(B,1)

            p_b = B(b,:);

            if ~isempty(B_ref) && size(B_ref,1) > 1
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
                p_b,x_line,hs_b,hf_ref_thick,rho_thick_const, ...
                hf_thick_high,rho_thick_high_manual);
        end

        rho_lo = prctile(rho_boot,2.5,2);
        rho_hi = prctile(rho_boot,97.5,2);

        hb = fill([x_line fliplr(x_line)], ...
            [rho_lo' fliplr(rho_hi')], ...
            [0.75 0.75 0.75], ...
            'FaceAlpha',0.25, ...
            'EdgeColor','none', ...
            'HandleVisibility','off');

        if i == 1
            h_band = hb;
        end

    end

    rho_line = hf_rho_model( ...
        p,x_line,hs_cases(i),hf_ref_thick,rho_thick_const, ...
        hf_thick_high,rho_thick_high_manual);

    h_lines(i) = plot(x_line,rho_line, ...
        'LineWidth',2.5, ...
        'Color',line_color, ...
        'DisplayName',hs_labels{i});
end

ylabel('Sea-ice density \rho (kg m^{-3})', ...
    'Interpreter','tex', ...
    'FontSize',11)

title(title_text, ...
    'Interpreter','tex', ...
    'FontWeight','normal', ...
    'FontSize',11)

if text_mode == "thresholds"
    hf_add_model_text_thresholds(gca,p,stats)
end

if doBootstrap && isgraphics(h_band)
    legend([h_lines; h_band; hFYI; hSYI], ...
        [hs_labels(:); {'95% CI'; 'FYI'; 'SYI & MYI'}], ...
        'Location','north', ...
        'NumColumns',3, ...
        'Box','off')
else
    legend([h_lines; hFYI; hSYI], ...
        [hs_labels(:); {'FYI'; 'SYI & MYI'}], ...
        'Location','north', ...
        'NumColumns',3, ...
        'Box','off')
end

xlim(x_lim)
ylim(rho_lim)
box on
grid off

end

function hf_add_model_text_thresholds(ax,p,stats)

hs_low = p(3);
hs_high = p(4);

wtxt = make_w_text(hs_low,hs_high);

txt = sprintf([ ...
    '$\\rho_{\\ell}(h_f),\\rho_h(h_f)\\;\\mathrm{fixed\\;from\\;reference\\;fit}$\n' ...
    '%s' ...
    '$R^2=%.2f,\\; RMSE=%.1f\\,\\mathrm{kg\\,m^{-3}},\\; N=%d$'], ...
    wtxt, ...
    stats.R2,stats.RMSE,stats.N);

text(ax,0.98,0.01,txt, ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom', ...
    'FontSize',10);

end

function c = hf_get_colormap_color(cmap,value,cmin,cmax)

value = min(max(value,cmin),cmax);

if cmax == cmin
    idx = 1;
else
    idx = 1 + round((value - cmin) ./ ...
        (cmax - cmin) .* (size(cmap,1) - 1));
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

function hf_print_fit(p,stats,hf_ref_thick,rho_thick_const,hf_thick_high,rho_thick_high_manual)

fprintf('low-snow rho at hf = 0          = %.2f kg/m3\n',p(1))
fprintf('low-snow constant rho           = %.2f kg/m3\n',p(2))
fprintf('low-snow hs threshold           = %.3f m\n',p(3))
fprintf('high-snow hs threshold          = %.3f m\n',p(4))
fprintf('critical hf low snow            = %.3f m\n',p(5))
fprintf('fixed rho_thick hf <= %.2f      = %.2f kg/m3\n',hf_ref_thick,rho_thick_const)
fprintf('manual rho_thick at hf = %.2f   = %.2f kg/m3\n',hf_thick_high,rho_thick_high_manual)
fprintf('RMSE = %.1f kg/m3\n',stats.RMSE)
fprintf('R2   = %.3f\n',stats.R2)

end
