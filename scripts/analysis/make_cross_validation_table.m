function make_cross_validation_table(repoRoot)
%
% MAKE_CROSS_VALIDATION_TABLE
%
% Performs leave-one-campaign-out cross-validation of the Arctic sea-ice
% density parameterization and generates Table S2.
%
% Workflow:
%   1. Load the merged density database.
%   2. Fit the temperature-thickness density model to all observations.
%   3. Refit the model repeatedly while withholding one campaign.
%   4. Predict densities for the withheld campaign.
%   5. Summarize campaign-specific and overall cross-validation skill.
%
% Input:
%   data/final/snow_model_matchups_with_models.mat
%
% Output:
%   results/TableS2_leave_one_campaign_out_cross_validation.csv
%

close all

finalDir = fullfile(repoRoot, 'data', 'final');
resultsDir = fullfile(repoRoot, 'results');

if ~exist(resultsDir, 'dir')
    mkdir(resultsDir)
end

mainFile = fullfile(finalDir, 'snow_model_matchups_with_models.mat');

load(mainFile, 'D')

densityMode = 1;

rho_cols = { ...
    "density_lab_kgm3", ...
    "density_insitu_connected_kgm3", ...
    "density_insitu_disconnected_kgm3"};

rho_names = {"lab","connected","disconnected"};

rho_col = rho_cols{densityMode};

fprintf('Using density mode: %s\n', rho_names{densityMode})

requiredVars = ["date","temperature_C","ice_thickness_m", ...
    "dataset","ice_age",rho_col];

missingVars = requiredVars(~ismember(requiredVars, ...
    string(D.Properties.VariableNames)));

if ~isempty(missingVars)
    error('D is missing required variables: %s', ...
        strjoin(cellstr(missingVars), ', '))
end

t = D.date;
T = D.temperature_C;
rho = D.(rho_col);
h = D.ice_thickness_m;
dataset_id = string(D.dataset);
ice_type = upper(strtrim(string(D.ice_age)));

if ~isdatetime(t)
    t = datetime(t);
end

valid = ...
    ~isnat(t) & ...
    isfinite(T) & ...
    isfinite(rho) & ...
    isfinite(h) & ...
    T < 0 & ...
    h > 0;

t = t(valid);
T = T(valid);
rho = rho(valid);
h = h(valid);
dataset_id = dataset_id(valid);
ice_type = ice_type(valid);

% Group individual dataset identifiers into campaign-level units used for leave-one-campaign-out cross-validation.
campaign = define_campaign_groups(dataset_id);

% Exclude unclassified datasets from campaign-level cross-validation.
campaign_list = unique(campaign);
campaign_list(campaign_list == "Other" | campaign_list == "") = [];

h_cold_peak = 1.10;
rho_cold_peak = 911;
h_cold_high = 3.50;
rho_cold_high = 908.5;

% Candidate parameter values explored during grid-search optimization.
rho_cold_0_grid = 905:1:912;
Tcrit_grid = -2.8:0.1:-2.1;
hcrit_warm_grid = 1.8:0.1:2.3;
rho_warm_plateau_grid = 895:1:905;

min_n_cold = 8;
min_n_warm = 8;

% Fit the full temperature-thickness density parameterization using all available observations.
fit = fit_fast_fixed_hcrit_model( ...
    rho, T, h, ...
    h_cold_peak, rho_cold_peak, h_cold_high, rho_cold_high, ...
    rho_cold_0_grid, Tcrit_grid, ...
    rho_warm_plateau_grid, hcrit_warm_grid, ...
    min_n_cold, min_n_warm);

p = [fit.Tcrit fit.rho_warm_plateau fit.warm_slope ...
     fit.hcrit_warm fit.rho_cold_0];

rhohat = predict_fast_fixed_hcrit_model( ...
    T, h, p, h_cold_peak, rho_cold_peak, ...
    h_cold_high, rho_cold_high);

metrics = regression_metrics(rho, rhohat, 5);

fprintf('\nFULL MODEL\n')
fprintf('-----------------------------------\n')
fprintf('Density mode = %s\n', rho_names{densityMode})
fprintf('Tcrit = %.2f deg C\n', p(1))
fprintf('rho_warm_plateau = %.1f kg m^-3\n', p(2))
fprintf('warm_slope = %.2f kg m^-3 m^-1\n', p(3))
fprintf('hcrit_warm = %.2f m\n', p(4))
fprintf('rho_cold_0 = %.1f kg m^-3\n', p(5))
fprintf('RMSE = %.2f kg m^-3\n', metrics.rmse)
fprintf('R2 = %.3f\n', metrics.r2)
fprintf('Adjusted R2 = %.3f\n', metrics.r2adj)
fprintf('N = %d\n', numel(rho))

rho_cv = nan(size(rho));
cv_summary = table();

% Leave-one-campaign-out cross-validation tests model transferability to entirely unseen observational campaigns.
for ic = 1:numel(campaign_list)

    campaign_name = campaign_list(ic);

    idx_test = campaign == campaign_name;
    idx_train = ~idx_test;

    if sum(idx_test) < 1 || sum(idx_train) < 10
        continue
    end

    try
        fit_cv = fit_fast_fixed_hcrit_model( ...
            rho(idx_train), T(idx_train), h(idx_train), ...
            h_cold_peak, rho_cold_peak, ...
            h_cold_high, rho_cold_high, ...
            rho_cold_0_grid, Tcrit_grid, ...
            rho_warm_plateau_grid, hcrit_warm_grid, ...
            min_n_cold, min_n_warm);
    catch ME
        warning('Skipping campaign %s: %s', campaign_name, ME.message)
        continue
    end

    p_cv = [fit_cv.Tcrit fit_cv.rho_warm_plateau ...
            fit_cv.warm_slope fit_cv.hcrit_warm ...
            fit_cv.rho_cold_0];

    rho_cv(idx_test) = predict_fast_fixed_hcrit_model( ...
        T(idx_test), h(idx_test), p_cv, ...
        h_cold_peak, rho_cold_peak, ...
        h_cold_high, rho_cold_high);

    resid = rho_cv(idx_test) - rho(idx_test);

    tmp = table();

    tmp.Campaign = campaign_name;
    tmp.N = sum(idx_test);
    tmp.MeanIceThickness_m = mean(h(idx_test), 'omitnan');
    tmp.RMSE_kgm3 = sqrt(mean(resid.^2, 'omitnan'));
    tmp.Bias_kgm3 = mean(resid, 'omitnan');

    cv_summary = [cv_summary; tmp];

end

idx_cv_valid = isfinite(rho_cv);

rmse_cv = sqrt(mean((rho(idx_cv_valid) - rho_cv(idx_cv_valid)).^2, ...
    'omitnan'));

bias_cv = mean(rho_cv(idx_cv_valid) - rho(idx_cv_valid), ...
    'omitnan');

rss_cv = sum((rho(idx_cv_valid) - rho_cv(idx_cv_valid)).^2, ...
    'omitnan');

sst_cv = sum((rho(idx_cv_valid) - mean(rho(idx_cv_valid), ...
    'omitnan')).^2, 'omitnan');

r2_cv = 1 - rss_cv / sst_cv;

fprintf('\nLEAVE-ONE-CAMPAIGN-OUT CROSS-VALIDATION\n')
fprintf('-----------------------------------\n')
fprintf('N_cv = %d\n', sum(idx_cv_valid))
fprintf('RMSE_CV = %.2f kg m^-3\n', rmse_cv)
fprintf('Bias_CV = %.2f kg m^-3\n', bias_cv)
fprintf('R2_CV = %.3f\n\n', r2_cv)

cv_summary = sortrows(cv_summary, 'Campaign');
cv_summary.MeanIceThickness_m = round(cv_summary.MeanIceThickness_m, 2);
cv_summary.RMSE_kgm3 = round(cv_summary.RMSE_kgm3, 1);
cv_summary.Bias_kgm3 = round(cv_summary.Bias_kgm3, 1);

disp(cv_summary)

outFile = fullfile(resultsDir, 'TableS2_leave_one_campaign_out_cross_validation.csv');
writetable(cv_summary, outFile)

fprintf('Saved cross-validation table to:\n%s\n', outFile)

end

% Grid-search optimization of the temperature-thickness density model.
function fit = fit_fast_fixed_hcrit_model(rho, T, h, ...
    h_cold_peak, rho_cold_peak, h_cold_high, rho_cold_high, ...
    rho_cold_0_grid, Tcrit_grid, rho_warm_plateau_grid, ...
    hcrit_warm_grid, min_n_cold, min_n_warm)

fit.rmse = inf;
fit.r2 = -inf;

sst = sum((rho - mean(rho, 'omitnan')).^2, 'omitnan');

for iT = 1:numel(Tcrit_grid)

    Tcrit = Tcrit_grid(iT);

    if sum(T <= Tcrit) < min_n_cold || sum(T > Tcrit) < min_n_warm
        continue
    end

    wT = ones(size(T));
    idxWarm = T > Tcrit;
    wT(idxWarm) = T(idxWarm) ./ Tcrit;
    wT = max(0, min(1, wT));

    for ih = 1:numel(hcrit_warm_grid)

        hcrit_warm = hcrit_warm_grid(ih);
        he_warm = min(h, hcrit_warm);
        Xh_warm = he_warm - hcrit_warm;

        for ic = 1:numel(rho_cold_0_grid)

            rho_cold_0 = rho_cold_0_grid(ic);

            rhoCold = cold_branch_density(h, rho_cold_0, ...
                h_cold_peak, rho_cold_peak, ...
                h_cold_high, rho_cold_high);

            for ir = 1:numel(rho_warm_plateau_grid)

                rho_warm_plateau = rho_warm_plateau_grid(ir);

                base = wT .* rhoCold + ...
                    (1 - wT) .* rho_warm_plateau;

                X = (1 - wT) .* Xh_warm;

                if all(abs(X) < eps)
                    continue
                end

                warm_slope = X \ (rho - base);

                if ~isfinite(warm_slope) || warm_slope < 0
                    continue
                end

                yhat = base + X .* warm_slope;

                rmse = sqrt(mean((rho - yhat).^2, 'omitnan'));
                r2 = 1 - sum((rho - yhat).^2, 'omitnan') / sst;

                if rmse < fit.rmse
                    fit.rmse = rmse;
                    fit.r2 = r2;
                    fit.Tcrit = Tcrit;
                    fit.rho_warm_plateau = rho_warm_plateau;
                    fit.warm_slope = warm_slope;
                    fit.hcrit_warm = hcrit_warm;
                    fit.rho_cold_0 = rho_cold_0;
                end

            end
        end
    end
end

if ~isfinite(fit.rmse)
    error('No valid model found. Widen grids or reduce min_n thresholds.')
end

end

% Evaluate the fitted density parameterization.
function rhohat = predict_fast_fixed_hcrit_model(T, h, p, ...
    h_cold_peak, rho_cold_peak, h_cold_high, rho_cold_high)

Tcrit = p(1);
rho_warm_plateau = p(2);
warm_slope = p(3);
hcrit_warm = p(4);
rho_cold_0 = p(5);

rhoCold = cold_branch_density(h, rho_cold_0, ...
    h_cold_peak, rho_cold_peak, ...
    h_cold_high, rho_cold_high);

he_warm = min(h, hcrit_warm);

rhoWarm = rho_warm_plateau + ...
    warm_slope .* (he_warm - hcrit_warm);

wT = ones(size(T));
idxWarm = T > Tcrit;
wT(idxWarm) = T(idxWarm) ./ Tcrit;
wT = max(0, min(1, wT));

rhohat = wT .* rhoCold + (1 - wT) .* rhoWarm;

end

% Piecewise-linear cold-regime density parameterization.
function rhoCold = cold_branch_density(h, rho_cold_0, ...
    h_cold_peak, rho_cold_peak, h_cold_high, rho_cold_high)

rhoCold = nan(size(h));

idx1 = h <= h_cold_peak;

rhoCold(idx1) = rho_cold_0 + ...
    (rho_cold_peak - rho_cold_0) .* h(idx1) ./ h_cold_peak;

idx2 = h > h_cold_peak;

rhoCold(idx2) = rho_cold_peak + ...
    (rho_cold_high - rho_cold_peak) .* ...
    (min(h(idx2), h_cold_high) - h_cold_peak) ./ ...
    (h_cold_high - h_cold_peak);

end

% Compute standard regression performance metrics.
function metrics = regression_metrics(y, yhat, k)

e = y - yhat;
n = numel(y);

sse = sum(e.^2, 'omitnan');
sst = sum((y - mean(y, 'omitnan')).^2, 'omitnan');

metrics.rmse = sqrt(mean(e.^2, 'omitnan'));
metrics.r2 = 1 - sse / sst;
metrics.r2adj = 1 - (1 - metrics.r2) * (n - 1) / max(n - k - 1, 1);

end

% Map dataset identifiers to campaign-level cross-validation groups.
function campaign_group = define_campaign_groups(dataset_id)

dataset_raw = strtrim(string(dataset_id));
dataset_id = lower(dataset_raw);

campaign_group = strings(size(dataset_id));
campaign_group(:) = "Other";

campaign_group(contains(dataset_id, "mosaic")) = "MOSAiC";
campaign_group(contains(dataset_id, "nansen")) = "Nansen Legacy";

campaign_group( ...
    contains(dataset_id, "gonorth") | ...
    contains(dataset_id, "go north") | ...
    contains(dataset_id, "go-north") | ...
    contains(dataset_id, "go_north") | ...
    contains(dataset_id, "gn")) = "GoNorth";

campaign_group( ...
    contains(dataset_id, "sudarco") | ...
    contains(dataset_id, "sudarko")) = "SUDARCO";

campaign_group(contains(dataset_id, "contrasts")) = "CONTRASTS";

campaign_group( ...
    contains(dataset_id, "svalbard") | ...
    contains(dataset_id, "kongsfjorden") | ...
    contains(dataset_id, "kongsfjord") | ...
    contains(dataset_id, "sva")) = "Svalbard";

unclassified = unique(dataset_raw(campaign_group == "Other"));

if ~isempty(unclassified)
    fprintf('\nUnclassified dataset values:\n')
    disp(unclassified)
end

end