function manuscript_statistics(repoRoot)

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

is_svalbard = is_svalbard(valid);

snow_obs = snow_obs(valid);
snow_era5 = snow_era5(valid);
snow_merra2 = snow_merra2(valid);

isFYI = strcmpi(age,'FYI');
isSYIMYI = strcmpi(age,'SYI') | strcmpi(age,'MYI');

%% Seasonal coverage

m = month(t);

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

idx = ~is_svalbard & ...
      isfinite(snow_obs) & ...
      isfinite(snow_era5);

x = snow_obs(idx);
y = snow_era5(idx);

N = numel(x);

RMSE = sqrt(mean((y - x).^2));
bias = mean(y - x);
bias_pct = 100 * bias / mean(x);

r = corr(x,y,'rows','complete');
R2 = r^2;

fprintf('ERA5    %5.2f  %5.2f   %6.3f   %+6.3f   %+7.1f   %d\n', ...
    r,R2,RMSE,bias,bias_pct,N)

idx = ~is_svalbard & ...
      isfinite(snow_obs) & ...
      isfinite(snow_merra2);

x = snow_obs(idx);
y = snow_merra2(idx);

N = numel(x);

RMSE = sqrt(mean((y - x).^2));
bias = mean(y - x);
bias_pct = 100 * bias / mean(x);

r = corr(x,y,'rows','complete');
R2 = r^2;

fprintf('MERRA2  %5.2f  %5.2f   %6.3f   %+6.3f   %+7.1f   %d\n', ...
    r,R2,RMSE,bias,bias_pct,N)

fprintf('\n')

%% Parameterization sample sizes

fprintf('Parameterization sample sizes\n')
fprintf('-----------------------------\n')

fprintf('Temperature-based parameterization N : %d\n', ...
    numel(rho))

fprintf('Snow-based parameterization N        : %d\n', ...
    sum(~is_svalbard))

fprintf('\n')

end