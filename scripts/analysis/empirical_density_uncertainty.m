function empirical_density_uncertainty(repoRoot)

if nargin < 1 || isempty(repoRoot)
    repoRoot = pwd;
end

close all

% make_empirical_density_uncertainty_new.m
%
% Estimate empirical sea-ice density uncertainty from observed residuals:
%
%   residual = rho_observed - rho_model(h_f, h_s,SM-LG)
%
% This version uses:
%   1. The updated data import used for Fig. 2:
%      data/final/snow_model_matchups_with_models.mat
%   2. The updated ERA5 SM-LG freeboard parametrization from Fig. 2d:
%
%      rho(h_f,h_s) = (1 - w) rho_l(h_f) + w rho_h(h_f)
%
%      rho_l = 865 + 129 h_f,        h_f < 0.28 m
%      rho_l = 900,                  h_f >= 0.28 m
%
%      rho_h = 911,                  h_f <= 0.20 m
%      rho_h = 911 - 20(h_f - 0.20), h_f > 0.20 m
%
%      w = h_s / 0.08,  0 < h_s < 0.08 m
%      w = 1,           h_s >= 0.08 m
%
% The table groups residuals by freeboard intervals and by the two
% snow regimes used for uncertainty reporting:
%
%   h_s < 0.08 m
%   h_s >= 0.08 m

finalDir = fullfile(repoRoot,'data','final');

inputFile = fullfile(finalDir,'snow_model_matchups_with_models.mat');

if ~exist(inputFile,'file')
    error('Input file not found: %s',inputFile)
end

load(inputFile,'D')

densityMode = 1; % 1 = lab, 2 = connected, 3 = disconnected

rho_cols = {"density_lab_kgm3", ...
            "density_insitu_connected_kgm3", ...
            "density_insitu_disconnected_kgm3"};

rho_names = {"lab","connected","disconnected"};
rho_col = rho_cols{densityMode};

fprintf('Using density mode: %s\n',rho_names{densityMode})

rho_w = 1024;

requiredVars = ["date","dataset","ice_thickness_m","hs_smlg",rho_col];
missingVars = requiredVars(~ismember(requiredVars,string(D.Properties.VariableNames)));

if ~isempty(missingVars)
    error('D is missing required variables: %s',strjoin(cellstr(missingVars),', '))
end

dataset_all = string(D.dataset);

if ismember("is_freeze",string(D.Properties.VariableNames))
    is_freeze = logical(D.is_freeze);
else
    is_freeze = false(height(D),1);
end

mask = ~(is_freeze | contains(lower(dataset_all),"svalbard"));

hi_all = D.ice_thickness_m(mask);
rho_obs_all = D.(rho_col)(mask);
date_all = D.date(mask);
hs_smlg_all = D.hs_smlg(mask);

[hf,hs_smlg,rho_obs] = prepare_freeboard_data( ...
    hi_all,rho_obs_all,hs_smlg_all,date_all,rho_w);

N = numel(rho_obs);

fprintf('Valid panel (d) data points: N = %d\n',N)

rho_model = rho_seaice_freeboard_snow_smlg_new(hf,hs_smlg);
residual = rho_obs - rho_model;

fb_edges = [0.00 0.10
            0.10 0.20
            0.20 0.30
            0.30 0.35];

fb_labels = [ ...
    "0.00-0.10"
    "0.10-0.20"
    "0.20-0.30"
    "0.30-0.35"];

hs_crit = 0.08;

snow_edges = [ ...
    -inf,    hs_crit
     hs_crit, inf];

snow_condition = [ ...
    "h_s < h_{s,crit}"
    "h_s >= h_{s,crit}"];

hs_ref_values = [0, hs_crit];

min_n_for_percentiles = 5;

rows = table;

for i = 1:size(fb_edges,1)

    hf_min = fb_edges(i,1);
    hf_max = fb_edges(i,2);
    hf_mid = mean([hf_min hf_max]);

    if i < size(fb_edges,1)
        in_fb = hf >= hf_min & hf < hf_max;
    else
        in_fb = hf >= hf_min & hf <= hf_max;
    end

    for j = 1:size(snow_edges,1)

        hs_min = snow_edges(j,1);
        hs_max = snow_edges(j,2);

        if isinf(hs_min)
            in_snow = hs_smlg < hs_max;
        elseif isinf(hs_max)
            in_snow = hs_smlg >= hs_min;
        else
            in_snow = hs_smlg >= hs_min & hs_smlg < hs_max;
        end

        idx = in_fb & in_snow & isfinite(residual);
        r = residual(idx);
        n_bin = numel(r);

        hs_ref = hs_ref_values(j);
        rho_ref = rho_seaice_freeboard_snow_smlg_new(hf_mid,hs_ref);

        Trow = table;

        Trow.fb_interval_m = fb_labels(i);
        Trow.hf_mid_m = hf_mid;
        Trow.snow_condition = snow_condition(j);
        Trow.hs_reference_m = hs_ref;
        Trow.n_observations = n_bin;
        Trow.rho_representative_kgm3 = rho_ref;

        if n_bin >= 1
            Trow.residual_bias_kgm3 = mean(r,'omitnan');
            Trow.residual_rmse_kgm3 = sqrt(mean(r.^2,'omitnan'));
            Trow.residual_std_kgm3 = std(r,'omitnan');
        else
            Trow.residual_bias_kgm3 = nan;
            Trow.residual_rmse_kgm3 = nan;
            Trow.residual_std_kgm3 = nan;
        end

        if n_bin >= min_n_for_percentiles
            res_p025 = prctile(r,2.5);
            res_p975 = prctile(r,97.5);
            note = "";
        else
            res_p025 = prctile(residual,2.5);
            res_p975 = prctile(residual,97.5);
            note = "pooled residual percentiles used";
        end

        Trow.residual_p025_kgm3 = res_p025;
        Trow.residual_p975_kgm3 = res_p975;
        Trow.rho_lower_kgm3 = rho_ref + res_p025;
        Trow.rho_upper_kgm3 = rho_ref + res_p975;
        Trow.note = note;

        rows = [rows; Trow];

    end
end

T = rows;

T.hf_mid_m = round(T.hf_mid_m,2);
T.hs_reference_m = round(T.hs_reference_m,2);
T.rho_representative_kgm3 = round(T.rho_representative_kgm3,1);
T.residual_bias_kgm3 = round(T.residual_bias_kgm3,1);
T.residual_rmse_kgm3 = round(T.residual_rmse_kgm3,1);
T.residual_std_kgm3 = round(T.residual_std_kgm3,1);
T.residual_p025_kgm3 = round(T.residual_p025_kgm3,1);
T.residual_p975_kgm3 = round(T.residual_p975_kgm3,1);
T.rho_lower_kgm3 = round(T.rho_lower_kgm3,1);
T.rho_upper_kgm3 = round(T.rho_upper_kgm3,1);

T_readable = table;

T_readable.Freeboard_interval = T.fb_interval_m;
T_readable.Freeboard_midpoint = string(compose('%.2f',T.hf_mid_m));
T_readable.Snow_regime = T.snow_condition;
T_readable.Snow_reference = string(compose('%.2f',T.hs_reference_m));
T_readable.N = string(T.n_observations);
T_readable.Density = string(compose('%.1f',T.rho_representative_kgm3));
T_readable.Density_p025 = string(compose('%.1f',T.rho_lower_kgm3));
T_readable.Density_p975 = string(compose('%.1f',T.rho_upper_kgm3));
T_readable.RMSE_residual = string(compose('%.1f',T.residual_rmse_kgm3));
T_readable.Bias_residual = string(compose('%.1f',T.residual_bias_kgm3));
T_readable.Note = string(T.note);
T_readable.Note(ismissing(T_readable.Note) | T_readable.Note == "NaN") = "";

nan_mask = T.n_observations == 0;
T_readable.RMSE_residual(nan_mask) = "N/A";
T_readable.Bias_residual(nan_mask) = "N/A";

units = table;
units.Freeboard_interval = "m";
units.Freeboard_midpoint = "m";
units.Snow_regime = "-";
units.Snow_reference = "m";
units.N = "-";
units.Density = "kg m^-3";
units.Density_p025 = "kg m^-3";
units.Density_p975 = "kg m^-3";
units.RMSE_residual = "kg m^-3";
units.Bias_residual = "kg m^-3";
units.Note = "-";

T_export = [units; T_readable];

outputCsv = fullfile(finalDir,'Empirical_density_uncertainty.csv');
writetable(T_export,outputCsv)

fprintf('Saved empirical density uncertainty table: %s\n',outputCsv)
fprintf('N = %d\n',N)
fprintf('Residual bias = %.1f kg m^-3\n',mean(residual,'omitnan'))
fprintf('Residual RMSE = %.1f kg m^-3\n',sqrt(mean(residual.^2,'omitnan')))
fprintf('Residual std  = %.1f kg m^-3\n',std(residual,'omitnan'))
fprintf('Number of uncertainty bins = %d\n',height(T))

%% Helpers

function [hf,hs,rho] = prepare_freeboard_data(hi,rho,hs,date,rho_w)

aug1_this_year = datetime(year(date),8,1);
aug1_next = aug1_this_year;
aug1_next(date > aug1_this_year) = datetime(year(date(date > aug1_this_year)) + 1,8,1);

rho_s = 0.35 .* days(aug1_next - date) + 239.78;

hf = hi - (rho.*hi + rho_s.*hs)./rho_w;
hf = max(hf,0);

hf = hf(:);
hs = hs(:);
rho = rho(:);

valid = isfinite(hf) & isfinite(hs) & isfinite(rho);

hf = hf(valid);
hs = hs(valid);
rho = rho(valid);

end

function rho_i = rho_seaice_freeboard_snow_smlg_new(hf,hs)

hf = double(hf);
hs = double(hs);

rho_low = 865 + 129 .* hf;
rho_low(hf >= 0.28) = 900;

rho_high = 911 .* ones(size(hf));
idx = hf > 0.20;
rho_high(idx) = 911 - 20 .* (hf(idx) - 0.20);

hs_crit = 0.08;
w = min(max(hs ./ hs_crit,0),1);

rho_i = (1 - w).*rho_low + w.*rho_high;

end

end
