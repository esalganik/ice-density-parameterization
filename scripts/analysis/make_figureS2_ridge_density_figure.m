function make_figureS2_ridge_density_figure(repoRoot)
%
% MAKE_SI_RIDGE_DENSITY_FIGURE
%
% Generates Supporting Figure S2 comparing directly measured MOSAiC ridge core densities with level-ice density observations.
%
% Input:
%   data/final/snow_model_matchups_with_models.mat
%   data/processed/Summary_Ridges.mat
%
% Output:
%   figures/FigS2.png
%
close all

finalDir = fullfile(repoRoot, 'data', 'final');
processedDir = fullfile(repoRoot, 'data', 'processed');
figureDir = fullfile(repoRoot, 'figures');
colormapDir = fullfile(repoRoot, 'data', 'colormaps');

if ~exist(figureDir, 'dir')
    mkdir(figureDir)
end

mainFile = fullfile(finalDir, 'snow_model_matchups_with_models.mat');
ridgeFile = fullfile(processedDir, 'Summary_Ridges.mat');

load(mainFile, 'D')
load(ridgeFile, 'Summary_Ridges')

S = load(fullfile(colormapDir, 'lipari.mat'));
lipari = S.lipari;

% Laboratory density is used for consistency with the main density parameterization figures.
rho_col = "density_lab_kgm3";

if ~isdatetime(D.date)
    D.date = datetime(D.date);
end

if ~isdatetime(Summary_Ridges.date)
    Summary_Ridges.date = datetime(Summary_Ridges.date);
end

T = D.temperature_C;
h = D.ice_thickness_m;
rho = D.(rho_col);
dataset_id = string(D.dataset);
ice_type = upper(string(D.ice_age));

is_svalbard = contains(lower(dataset_id), "svalbard");

valid = isfinite(T) & isfinite(h) & isfinite(rho) & ...
    T < 0 & h > 0;

T = T(valid);
h = h(valid);
rho = rho(valid);
ice_type = ice_type(valid);
is_svalbard = is_svalbard(valid);

% Separate level-ice observations into FYI, older ice, and Svalbard groups for background context behind the ridge-core measurements.
is_fyi = ice_type == "FYI" & ~is_svalbard;
is_old = ice_type ~= "FYI" & ~is_svalbard;
is_sva = is_svalbard;

T_ridge = Summary_Ridges.temperature_C;
h_ridge = Summary_Ridges.ice_thickness_m;
rho_ridge = Summary_Ridges.(rho_col);

validRidge = isfinite(T_ridge) & isfinite(h_ridge) & ...
    isfinite(rho_ridge) & T_ridge < 0 & h_ridge > 0;

T_ridge = T_ridge(validRidge);
h_ridge = h_ridge(validRidge);
rho_ridge = rho_ridge(validRidge);

fs_ax = 10.5;
fs_txt = 9.0;
edgeCol = [0.75 0.75 0.75];
fillCol = [0.87 0.87 0.87];

ms = 32;
ms_ridge = 90;
lw_open = 0.9;
lw_ridge = 1.3;

hmin = 0;
hmax = 6;

Tmin = floor(min(T_ridge, [], 'omitnan') * 2) / 2;
Tmax = ceil(max(T_ridge, [], 'omitnan') * 2) / 2;

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1.0 1.0 10.8 4.2]);

tl = tiledlayout(fig, 1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax1 = nexttile(tl, 1);
hold(ax1, 'on')

scatter(ax1, T(is_fyi), rho(is_fyi), ms, ...
    'o', ...
    'MarkerFaceColor', fillCol, ...
    'MarkerEdgeColor', edgeCol, ...
    'LineWidth', 0.5)

scatter(ax1, T(is_old), rho(is_old), ms, ...
    'o', ...
    'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', edgeCol, ...
    'LineWidth', lw_open)

scatter(ax1, T(is_sva), rho(is_sva), ms, ...
    '>', ...
    'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', edgeCol, ...
    'LineWidth', lw_open)

scatter(ax1, T_ridge, rho_ridge, ms_ridge, h_ridge, ...
    'd', 'filled', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', lw_ridge)

colormap(ax1, lipari)
clim(ax1, [hmin hmax])

cb1 = colorbar(ax1);
ylabel(cb1, 'Ridge thickness (m)', 'FontSize', fs_ax)

xlabel(ax1, 'Sea-ice temperature{\it T_i} (°C)', 'Interpreter', 'tex')
ylabel(ax1, 'Sea-ice density \rho (kg m^{-3})', 'Interpreter', 'tex')
xlim(ax1, [-10 0])
ylim(ax1, [850 920])
box(ax1, 'on')
set(ax1, 'FontSize', fs_ax, 'FontWeight', 'normal')
text(ax1, -0.145, 1.03, '(a)', 'Units', 'normalized', ...
    'FontSize', fs_ax, 'FontWeight', 'normal')

ax2 = nexttile(tl, 2);
hold(ax2, 'on')

scatter(ax2, h(is_fyi), rho(is_fyi), ms, ...
    'o', ...
    'MarkerFaceColor', fillCol, ...
    'MarkerEdgeColor', edgeCol, ...
    'LineWidth', 0.5)

scatter(ax2, h(is_old), rho(is_old), ms, ...
    'o', ...
    'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', edgeCol, ...
    'LineWidth', lw_open)

scatter(ax2, h(is_sva), rho(is_sva), ms, ...
    '>', ...
    'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', edgeCol, ...
    'LineWidth', lw_open)

scatter(ax2, h_ridge, rho_ridge, ms_ridge, T_ridge, ...
    'd', 'filled', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', lw_ridge)

colormap(ax2, lipari)
clim(ax2, [Tmin Tmax])

cb2 = colorbar(ax2);
ylabel(cb2, 'Ridge temperature (°C)', 'FontSize', fs_ax)

xlabel(ax2, 'Sea-ice thickness{\it h_i} (m)', 'Interpreter', 'tex')
ylabel(ax2, 'Sea-ice density \rho (kg m^{-3})', 'Interpreter', 'tex')
xlim(ax2, [0 6])
ylim(ax2, [850 920])
box(ax2, 'on')
set(ax2, 'FontSize', fs_ax, 'FontWeight', 'normal')
text(ax2, -0.145, 1.03, '(b)', 'Units', 'normalized', ...
    'FontSize', fs_ax, 'FontWeight', 'normal')

p_fyi = plot(ax2, nan, nan, 'o', ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', fillCol, ...
    'MarkerEdgeColor', edgeCol, ...
    'LineStyle', 'none');

p_old = plot(ax2, nan, nan, 'o', ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', edgeCol, ...
    'LineWidth', lw_open, ...
    'LineStyle', 'none');

p_sva = plot(ax2, nan, nan, '>', ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', edgeCol, ...
    'LineWidth', lw_open, ...
    'LineStyle', 'none');

p_ridge = plot(ax2, nan, nan, 'd', ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'LineStyle', 'none');

legend(ax2, [p_fyi p_old p_sva p_ridge], ...
    {'FYI', 'SYI & MYI', 'Svalbard', 'Ridges'}, ...
    'Location', 'southeast', ...
    'Box', 'off', ...
    'FontSize', fs_txt, ...
    'Interpreter', 'tex')

set(findall(fig, 'Type', 'axes'), 'Toolbar', [])

outputFigure = fullfile(figureDir, 'FigS2.png');
exportgraphics(fig, outputFigure, 'Resolution', 300)

close(fig)

fprintf('Generated SI ridge density figure.\n')
fprintf('Main observations: %d\n', numel(rho))
fprintf('Ridge observations: %d\n', numel(rho_ridge))
fprintf('Saved figure to:\n%s\n', outputFigure)

end