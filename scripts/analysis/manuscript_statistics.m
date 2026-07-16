function manuscript_statistics(repoRoot)
%
% MANUSCRIPT_STATISTICS
%
% Generates summary statistics reported in the manuscript and supporting
% information.
%
% Outputs include:
%   - Dataset coverage and sample sizes
%   - Seasonal and regional observation statistics
%   - Ice-age-specific density summaries
%   - SnowModel-LG snow-thickness evaluation metrics
%   - Parameterization sample sizes
%
% Input:
%   data/final/snow_model_matchups_with_models.mat
%

close all

finalDir = fullfile(repoRoot, ...
    'data', ...
    'final');

inputFile = fullfile(finalDir, ...
    'snow_model_matchups_with_models.mat');

load(inputFile,'D')

densityMode = 1; % 1 = lab, 2 = connected, 3 = disconnected

rho_cols = {"density_lab_kgm3", ...
            "density_insitu_connected_kgm3", ...
            "density_insitu_disconnected_kgm3"};

rho_col = rho_cols{densityMode};

requiredVars = ["date", ...
                "dataset", ...
                "ice_age", ...
                "temperature_C", ...
                "ice_thickness_m", ...
                "hs", ...
                "hs_smlg", ...
                "hs_merra", ...
                rho_col];

missingVars = requiredVars( ...
    ~ismember(requiredVars, ...
    string(D.Properties.VariableNames)));

if ~isempty(missingVars)

    error('D is missing required variables: %s', ...
        strjoin(cellstr(missingVars),', '))

end

if ~isdatetime(D.date)
    D.date = datetime(D.date);
end

t = D.date;
T = D.temperature_C;
rho = D.(rho_col);
h = D.ice_thickness_m;
age = string(D.ice_age);
dataset_id = string(D.dataset);

snow_obs = D.hs;
snow_era5 = D.hs_smlg;
snow_merra2 = D.hs_merra;

is_svalbard = contains(lower(dataset_id),"svalbard");

% Retain only physically valid sea-ice observations used throughout the parameterization analysis.
valid = ~isnat(t) & ...
        isfinite(T) & ...
        isfinite(rho) & ...
        isfinite(h) & ...
        T < 0 & ...
        h > 0;

t = t(valid);
T = T(valid);
rho = rho(valid);
h = h(valid);
age = age(valid);
dataset_id = dataset_id(valid);

is_svalbard = is_svalbard(valid);

snow_obs = snow_obs(valid);
snow_era5 = snow_era5(valid);
snow_merra2 = snow_merra2(valid);

isFYI = strcmpi(age,'FYI');
isSYIMYI = strcmpi(age,'SYI') | strcmpi(age,'MYI');

%% Dataset overview

fprintf('\n')
fprintf('Dataset overview\n')
fprintf('----------------\n')
fprintf(['%-14s %-11s %-22s %-20s %-16s %-18s %4s\n'], ...
    'Dataset','Years','Months','Region','Ice type','Ice thickness (m)','N')

datasetOrder = ["N-ICE2015", ...
                "MOSAiC", ...
                "GoNorth", ...
                "SUDARCO", ...
                "CONTRASTS", ...
                "Nansen Legacy", ...
                "Svalbard"];

for i = 1:numel(datasetOrder)

    ds = datasetOrder(i);
    switch ds
    case "N-ICE2015"
        idx = contains(dataset_id,"N-ICE","IgnoreCase",true);
    case "MOSAiC"
        idx = contains(dataset_id,"MOSAiC","IgnoreCase",true);
    case "GoNorth"
        idx = contains(dataset_id,"GoNorth","IgnoreCase",true);
    case "SUDARCO"
        idx = contains(dataset_id,"SUDARCO","IgnoreCase",true);
    case "CONTRASTS"
        idx = contains(dataset_id,"CONTRASTS","IgnoreCase",true);
    case "Nansen Legacy"
        idx = contains(dataset_id,"Nansen","IgnoreCase",true);
    case "Svalbard"
        idx = contains(dataset_id,"Svalbard","IgnoreCase",true);
    end

    if ~any(idx)
        continue
    end

    yearsStr = format_years_for_table(year(t(idx)));
    monthsStr = format_months_for_table(month(t(idx)));
    iceTypeStr = strjoin(unique(age(idx))',', ');

    if ds == "Svalbard"
        regionStr = "Kongsfjorden";
    elseif ds == "Nansen Legacy"
        regionStr = "AO, Barents Sea";
    else
        regionStr = "Arctic Ocean (AO)";
    end

    fprintf('%-14s %-11s %-22s %-20s %-16s %.1f--%.1f          %4d\n', ...
        ds, ...
        yearsStr, ...
        monthsStr, ...
        regionStr, ...
        iceTypeStr, ...
        min(h(idx)), ...
        max(h(idx)), ...
        sum(idx))

end

fprintf('\n')

%% Seasonal coverage

m = month(t);

% Meteorological seasons.
season = strings(size(m));

season(ismember(m,[12 1 2])) = "Winter";
season(ismember(m,[3 4 5])) = "Spring";
season(ismember(m,[6 7 8])) = "Summer";
season(ismember(m,[9 10 11])) = "Autumn";

fprintf('\n')
fprintf('Seasonal coverage\n')
fprintf('-----------------\n')
fprintf('            Winter  Spring  Summer  Autumn   Total\n')

fprintf('Arctic      %6d  %6d  %6d  %6d  %6d\n', ...
    sum(~is_svalbard & season=="Winter"), ...
    sum(~is_svalbard & season=="Spring"), ...
    sum(~is_svalbard & season=="Summer"), ...
    sum(~is_svalbard & season=="Autumn"), ...
    sum(~is_svalbard))

fprintf('Svalbard    %6d  %6d  %6d  %6d  %6d\n', ...
    sum(is_svalbard & season=="Winter"), ...
    sum(is_svalbard & season=="Spring"), ...
    sum(is_svalbard & season=="Summer"), ...
    sum(is_svalbard & season=="Autumn"), ...
    sum(is_svalbard))

fprintf('\n')

%% Dataset characteristics

fprintf('Dataset characteristics\n')
fprintf('-----------------------\n')

fprintf('Total valid observations     : %d\n',numel(rho))
fprintf('Arctic observations          : %d\n',sum(~is_svalbard))
fprintf('Svalbard observations        : %d\n',sum(is_svalbard))

fprintf('\n')

fprintf('FYI                          : %d\n',sum(isFYI))
fprintf('SYI/MYI                      : %d\n',sum(isSYIMYI))

fprintf('\n')

fprintf('Ice thickness (m)\n')
fprintf('  Range       : %.2f - %.2f\n', ...
    min(h),max(h))
fprintf('  Mean +/- SD : %.2f +/- %.2f\n', ...
    mean(h),std(h))

fprintf('\n')

fprintf('Temperature (deg C)\n')
fprintf('  Range       : %.1f - %.1f\n', ...
    min(T),max(T))
fprintf('  Mean +/- SD : %.1f +/- %.1f\n', ...
    mean(T),std(T))

fprintf('\n')

fprintf('Density (kg m^-3)\n')
fprintf('  Range       : %.0f - %.0f\n', ...
    min(rho),max(rho))
fprintf('  Mean +/- SD : %.0f +/- %.0f\n', ...
    mean(rho),std(rho))

fprintf('\n')

fprintf('Observed snow thickness (m)\n')
fprintf('  Range       : %.2f - %.2f\n', ...
    min(snow_obs,[],'omitnan'), ...
    max(snow_obs,[],'omitnan'))
fprintf('  Mean +/- SD : %.2f +/- %.2f\n', ...
    mean(snow_obs,'omitnan'), ...
    std(snow_obs,'omitnan'))

fprintf('\n')

%% Ice-age summary

fprintf('Ice-age density summary\n')
fprintf('-----------------------\n')

fprintf('FYI      : %.0f +/- %.0f kg m^-3 (N = %d)\n', ...
    mean(rho(isFYI),'omitnan'), ...
    std(rho(isFYI),'omitnan'), ...
    sum(isFYI))

fprintf('SYI/MYI  : %.0f +/- %.0f kg m^-3 (N = %d)\n', ...
    mean(rho(isSYIMYI),'omitnan'), ...
    std(rho(isSYIMYI),'omitnan'), ...
    sum(isSYIMYI))

fprintf('\n')

%% Regional summary

fprintf('Regional density summary\n')
fprintf('------------------------\n')

fprintf('Arctic   : %.0f +/- %.0f kg m^-3 (N = %d)\n', ...
    mean(rho(~is_svalbard),'omitnan'), ...
    std(rho(~is_svalbard),'omitnan'), ...
    sum(~is_svalbard))

fprintf('Svalbard : %.0f +/- %.0f kg m^-3 (N = %d)\n', ...
    mean(rho(is_svalbard),'omitnan'), ...
    std(rho(is_svalbard),'omitnan'), ...
    sum(is_svalbard))

fprintf('\n')

%% SnowModel-LG evaluation

fprintf('Snow thickness models vs observations (non-Svalbard)\n')
fprintf('----------------------------------------------------\n')
fprintf('Model     r      R2     RMSE_m   Bias_m   Bias_pct    N\n')

model_names = {'ERA5','MERRA2'};
model_snow = {snow_era5, snow_merra2};

for k = 1:numel(model_names)

    % Restrict evaluation to Arctic pack-ice observations because SnowModel-LG is not intended to represent Svalbard fjord conditions.
    idx = ~is_svalbard & ...
          isfinite(snow_obs) & ...
          isfinite(model_snow{k});

    x = snow_obs(idx);
    y = model_snow{k}(idx);

    N = numel(x);
    RMSE = sqrt(mean((y - x).^2));
    bias = mean(y - x);
    bias_pct = 100 * bias / mean(x);

    r = corr(x,y,'rows','complete');
    R2 = r^2;

    fprintf('%-7s %5.2f  %5.2f   %6.3f   %+6.3f   %+7.1f   %d\n', ...
        model_names{k},r,R2,RMSE,bias,bias_pct,N)

end

fprintf('\n')

%% SnowModel-LG seasonal snow-bias evaluation

fprintf('SnowModel-LG snow bias by season group (non-Svalbard)\n')
fprintf('-----------------------------------------------------\n')
fprintf('Group                 Model     MeanObs  MeanModel  Bias_m  Bias_cm  Bias_pct  RMSE_m   R2     N\n')

groups = { ...
    true(size(season)), ...
    season ~= "Summer", ...
    season == "Winter", ...
    season == "Spring", ...
    season == "Summer", ...
    season == "Autumn"};

group_names = { ...
    'All seasons', ...
    'Non-summer', ...
    'Winter', ...
    'Spring', ...
    'Summer', ...
    'Autumn'};

for g = 1:numel(groups)

    for k = 1:numel(model_names)

        idx = ~is_svalbard & ...
              groups{g} & ...
              isfinite(snow_obs) & ...
              isfinite(model_snow{k});

        x = snow_obs(idx);
        y = model_snow{k}(idx);

        N = numel(x);

        if N > 1
            mean_obs = mean(x);
            mean_model = mean(y);
            bias = mean(y - x);
            bias_cm = 100 * bias;
            bias_pct = 100 * bias / mean_obs;
            RMSE = sqrt(mean((y - x).^2));
            r = corr(x,y,'rows','complete');
            R2 = r^2;
        else
            mean_obs = NaN;
            mean_model = NaN;
            bias = NaN;
            bias_cm = NaN;
            bias_pct = NaN;
            RMSE = NaN;
            R2 = NaN;
        end

        fprintf('%-20s  %-7s  %7.3f  %9.3f  %+6.3f  %+7.1f  %+8.1f  %6.3f  %5.2f  %4d\n', ...
            group_names{g}, ...
            model_names{k}, ...
            mean_obs, ...
            mean_model, ...
            bias, ...
            bias_cm, ...
            bias_pct, ...
            RMSE, ...
            R2, ...
            N)

    end

end

fprintf('\n')

%% Empirical density uncertainty

fprintf('Empirical density uncertainty\n')
fprintf('-----------------------------\n')

empirical_density_uncertainty(repoRoot)

fprintf('\n')

%% Parameterization sample sizes

fprintf('Parameterization sample sizes\n')
fprintf('-----------------------------\n')

fprintf('Temperature-based parameterization N : %d\n', ...
    numel(rho))

fprintf('Snow-based parameterization N        : %d\n', ...
    sum(~is_svalbard))

fprintf('\n')

% Format year ranges for manuscript summary tables.
function yearsStr = format_years_for_table(years)

years = unique(sort(years(:)'));

if isscalar(years)
    yearsStr = string(years);
else
    yearsStr = string(years(1)) + "--" + string(years(end));
end

end

% Format month ranges for manuscript summary tables.
function monthsStr = format_months_for_table(months)

monthNames = ["Jan","Feb","Mar","Apr","May","Jun", ...
              "Jul","Aug","Sep","Oct","Nov","Dec"];

months = unique(sort(months(:)'));

if isscalar(months)
    monthsStr = monthNames(months);
elseif all(diff(months) == 1)
    monthsStr = monthNames(months(1)) + "--" + monthNames(months(end));
else
    monthsStr = strjoin(monthNames(months),", ");
end

end

end