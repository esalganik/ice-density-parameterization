function make_figureS6_snow_model_evaluation(repoRoot)
%
% MAKE_FIGURES6_SNOW_MODEL_EVALUATION
%
% Generates Supporting Figure S6 comparing observed snow thickness with
% SnowModel-LG snow thickness forced by ERA5 and MERRA-2.
%
% The figure includes all non-Svalbard observations and colors points by
% season. Skill statistics R2 and RMSE are computed using all seasons,
% whereas the bias shown in the figure is computed for the non-summer season
% to match the manuscript text.
%
% Input:
%   data/final/snow_model_matchups_with_models.mat
%
% Output:
%   figures/FigS6.png
%

close all

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

%% Settings

densityMode = 1; % 1 = lab, 2 = connected, 3 = disconnected

%% Paths

finalDir = fullfile(repoRoot, ...
    'data', ...
    'final');

figureDir = fullfile(repoRoot, ...
    'figures');

if ~exist(figureDir,'dir')
    mkdir(figureDir)
end

inputFile = fullfile(finalDir, ...
    'snow_model_matchups_with_models.mat');

outputFigure = fullfile(figureDir, ...
    'FigS6.png');

%% Load data

load(inputFile,'D')

rho_cols = {"density_lab_kgm3", ...
            "density_insitu_connected_kgm3", ...
            "density_insitu_disconnected_kgm3"};

rho_names = {"lab","connected","disconnected"};
rho_col = rho_cols{densityMode};

requiredVars = ["date", ...
                "dataset", ...
                "temperature_C", ...
                "ice_thickness_m", ...
                "hs", ...
                "hs_smlg", ...
                "hs_merra", ...
                rho_col];

missingVars = requiredVars(~ismember(requiredVars,string(D.Properties.VariableNames)));

if ~isempty(missingVars)
    error('D is missing required variables: %s',strjoin(cellstr(missingVars),', '))
end

if ~isdatetime(D.date)
    D.date = datetime(D.date);
end

t = D.date;
T = D.temperature_C;
rho = D.(rho_col);
h = D.ice_thickness_m;
dataset_id = string(D.dataset);

snow_obs = D.hs;
snow_era5 = D.hs_smlg;
snow_merra2 = D.hs_merra;

is_svalbard = contains(lower(dataset_id),"svalbard");

%% Basic validity filter

valid = ~isnat(t) & ...
        isfinite(T) & ...
        isfinite(rho) & ...
        isfinite(h) & ...
        T < 0 & ...
        h > 0;

t = t(valid);
snow_obs = snow_obs(valid);
snow_era5 = snow_era5(valid);
snow_merra2 = snow_merra2(valid);
is_svalbard = is_svalbard(valid);

%% Season assignment

m = month(t);

season = strings(size(m));
season(ismember(m,[12 1 2])) = "Winter";
season(ismember(m,[3 4 5])) = "Spring";
season(ismember(m,[6 7 8])) = "Summer";
season(ismember(m,[9 10 11])) = "Autumn";

subsetLabel = 'all seasons, non-Svalbard';

analysisMask = ~is_svalbard;
nonSummerMask = analysisMask & season ~= "Summer";

%% Matchup statistics

% All-season statistics for R2, RMSE, and N.
idx_era5 = analysisMask & isfinite(snow_obs) & isfinite(snow_era5);
idx_merra = analysisMask & isfinite(snow_obs) & isfinite(snow_merra2);

obs_era5 = snow_obs(idx_era5);
mod_era5 = snow_era5(idx_era5);
season_era5 = season(idx_era5);

obs_merra = snow_obs(idx_merra);
mod_merra = snow_merra2(idx_merra);
season_merra = season(idx_merra);

stats_era5 = snow_metrics(obs_era5,mod_era5);
stats_merra = snow_metrics(obs_merra,mod_merra);

% Non-summer bias only, matching the manuscript wording.
idx_era5_ns = nonSummerMask & isfinite(snow_obs) & isfinite(snow_era5);
idx_merra_ns = nonSummerMask & isfinite(snow_obs) & isfinite(snow_merra2);

stats_era5_ns = snow_metrics(snow_obs(idx_era5_ns),snow_era5(idx_era5_ns));
stats_merra_ns = snow_metrics(snow_obs(idx_merra_ns),snow_merra2(idx_merra_ns));

fprintf('\n')
fprintf('Snow thickness models vs observations: all-season skill, non-summer bias\n')
fprintf('-----------------------------------------------------------------------\n')
fprintf('Density mode: %s\n',rho_names{densityMode})
fprintf('Model      R2_all  RMSE_all_m  N_all   Bias_nonSummer_m  Bias_cm  Bias_pct  N_nonSummer\n')

fprintf('%-8s %6.2f    %7.3f   %5d       %+7.3f     %+6.1f    %+7.1f     %5d\n', ...
    'ERA5', ...
    stats_era5.R2, ...
    stats_era5.RMSE, ...
    stats_era5.N, ...
    stats_era5_ns.bias, ...
    stats_era5_ns.bias_cm, ...
    stats_era5_ns.bias_pct, ...
    stats_era5_ns.N)

fprintf('%-8s %6.2f    %7.3f   %5d       %+7.3f     %+6.1f    %+7.1f     %5d\n', ...
    'MERRA-2', ...
    stats_merra.R2, ...
    stats_merra.RMSE, ...
    stats_merra.N, ...
    stats_merra_ns.bias, ...
    stats_merra_ns.bias_cm, ...
    stats_merra_ns.bias_pct, ...
    stats_merra_ns.N)

fprintf('\n')

%% Figure

fig = figure('Color','w');
set(fig,'Units','inches','Position',[1 1 12.0 5.2])

tiledlayout(1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Common axis limits for both panels.
maxSnow = max([obs_era5(:); mod_era5(:); obs_merra(:); mod_merra(:)],[],'omitnan');
snowLim = [0 ceil(maxSnow*10)/10];

if snowLim(2) < 0.5
    snowLim(2) = 0.5;
end

% Four-season colors.
seasonColors.Winter = [0.10 0.25 0.45];
seasonColors.Spring = [0.35 0.60 0.75];
seasonColors.Summer = [0.90 0.75 0.35];
seasonColors.Autumn = [0.80 0.45 0.25];

% Panel a: ERA5
ax1 = nexttile;
plot_snow_scatter(ax1, ...
    obs_era5,mod_era5,season_era5, ...
    stats_era5,stats_era5_ns, ...
    'ERA5 SM-LG', ...
    snowLim,subsetLabel,seasonColors,true)

text(ax1,-0.14,1.035,'(a)', ...
    'Units','normalized', ...
    'FontSize',11, ...
    'FontWeight','normal')

% Panel b: MERRA-2
ax2 = nexttile;
plot_snow_scatter(ax2, ...
    obs_merra,mod_merra,season_merra, ...
    stats_merra,stats_merra_ns, ...
    'MERRA-2 SM-LG', ...
    snowLim,subsetLabel,seasonColors,true)

text(ax2,-0.14,1.035,'(b)', ...
    'Units','normalized', ...
    'FontSize',11, ...
    'FontWeight','normal')

set(findall(fig,'Type','axes'), ...
    'FontSize',11, ...
    'LineWidth',1, ...
    'Box','on', ...
    'TickDir','in')

set(findall(fig,'Type','axes'),'Toolbar',[])

exportgraphics(fig,outputFigure,'Resolution',300)

close(fig)

fprintf('Generated Figure S6 using %s density.\n', ...
    rho_names{densityMode})

fprintf('Saved figure to:\n%s\n', ...
    outputFigure)

%% Helper functions

function stats = snow_metrics(obs,model)

obs = obs(:);
model = model(:);

idx = isfinite(obs) & isfinite(model);
obs = obs(idx);
model = model(idx);

stats.N = numel(obs);

if stats.N > 1

    stats.mean_obs = mean(obs,'omitnan');
    stats.mean_model = mean(model,'omitnan');
    stats.bias = mean(model - obs,'omitnan');
    stats.bias_cm = 100 .* stats.bias;
    stats.bias_pct = 100 .* stats.bias ./ stats.mean_obs;
    stats.RMSE = sqrt(mean((model - obs).^2,'omitnan'));

    % Same R2 definition as manuscript_statistics:
    % squared Pearson correlation.
    stats.r = corr(obs,model,'rows','complete');
    stats.R2 = stats.r.^2;

else

    stats.mean_obs = NaN;
    stats.mean_model = NaN;
    stats.bias = NaN;
    stats.bias_cm = NaN;
    stats.bias_pct = NaN;
    stats.RMSE = NaN;
    stats.r = NaN;
    stats.R2 = NaN;

end

end

function plot_snow_scatter(ax,obs,model,season,statsAll,statsNonSummer, ...
    titleText,snowLim,subsetLabel,seasonColors,showLegend)

hold(ax,'on')
box(ax,'on')

% 1:1 line behind points.
plot(ax,snowLim,snowLim, ...
    '--', ...
    'Color',[0.15 0.15 0.15], ...
    'LineWidth',1.2, ...
    'HandleVisibility','off')

seasonOrder = ["Winter","Spring","Summer","Autumn"];
seasonLabels = ["Winter","Spring","Summer","Autumn"];

hSeason = gobjects(numel(seasonOrder),1);

for i = 1:numel(seasonOrder)

    idx = season == seasonOrder(i);

    if seasonOrder(i) == "Winter"
        c = seasonColors.Winter;
    elseif seasonOrder(i) == "Spring"
        c = seasonColors.Spring;
    elseif seasonOrder(i) == "Summer"
        c = seasonColors.Summer;
    else
        c = seasonColors.Autumn;
    end

    hSeason(i) = scatter(ax,obs(idx),model(idx),38, ...
        'filled', ...
        'MarkerFaceColor',c, ...
        'MarkerFaceAlpha',0.82, ...
        'MarkerEdgeColor',[0.35 0.35 0.35], ...
        'LineWidth',0.6, ...
        'DisplayName',seasonLabels(i));

end

xlim(ax,snowLim)
ylim(ax,snowLim)
axis(ax,'square')

xlabel(ax,'Observed snow thickness (m)','FontSize',11)
ylabel(ax,'SM-LG snow thickness (m)','FontSize',11)
title(ax,titleText,'FontWeight','normal','FontSize',11)

txt = sprintf([ ...
    '%s\n' ...
    'N = %d\n' ...
    'R^2 = %.2f\n' ...
    'RMSE = %.3f m\n' ...
    'Non-summer bias = %+0.3f m (%+.0f%%)'], ...
    subsetLabel, ...
    statsAll.N, ...
    statsAll.R2, ...
    statsAll.RMSE, ...
    statsNonSummer.bias, ...
    statsNonSummer.bias_pct);

text(ax,0.96,0.02,txt, ...
    'Units','normalized', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','right', ...
    'FontSize',9.3)

if showLegend
    legend(ax,hSeason,seasonLabels, ...
        'Location','east', ...
        'Box','off', ...
        'FontSize',9.5)
end

end

end