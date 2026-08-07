%
% Reproduction checks for CryoSat-2 values reported in the manuscript.
%
% Produces:
%   1. Mean pan-Arctic sea-ice volume for February and July using the
%      A2010, J2022, and This Study density parameterizations.
%   2. Relative differences in sea-ice volume between parameterizations.
%   3. Climatological sea-ice thickness differences used to check the
%      regional CryoSat-2 statements in the manuscript.
%
% This script is intended only as a numerical manuscript check. It does
% not generate Figures 3 or S5.
%
% Required input files:
%   data/cs2/uit_sit_v2northpolarstereo_80km_feb_2011-2023_v3p0.nc
%   data/cs2/uit_sit_v2northpolarstereo_80km_jul_2011-2023_v3p0.nc
%

clear
clc

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));

cs2Dir = fullfile(repoRoot,'data','cs2');

filename_feb = fullfile(cs2Dir, ...
    'uit_sit_v2northpolarstereo_80km_feb_2011-2023_v3p0.nc');

filename_jul = fullfile(cs2Dir, ...
    'uit_sit_v2northpolarstereo_80km_jul_2011-2023_v3p0.nc');

if ~isfile(filename_feb)
    error('Missing CryoSat-2 input file: %s',filename_feb)
end

if ~isfile(filename_jul)
    error('Missing CryoSat-2 input file: %s',filename_jul)
end

% Read only the variables required for the manuscript checks.
%
% The NetCDF variables retain the internal J2021 and S2026 names used in
% the CryoSat-2 files. In the manuscript these parameterizations are
% referred to as J2022 and This Study, respectively.

grid_lat = ncread(filename_feb,'latitude');

sic_feb = ncread(filename_feb,'sic');
sic_jul = ncread(filename_jul,'sic');

sit_A2010_feb = ncread(filename_feb,'sit_A2010');
sit_J2022_feb = ncread(filename_feb,'sit_J2021');
sit_Study_feb = ncread(filename_feb,'sit_S2026');

sit_A2010_jul = ncread(filename_jul,'sit_A2010');
sit_J2022_jul = ncread(filename_jul,'sit_J2021');
sit_Study_jul = ncread(filename_jul,'sit_S2026');

siv_A2010_feb = ncread(filename_feb,'siv_A2010');
siv_J2022_feb = ncread(filename_feb,'siv_J2021');
siv_Study_feb = ncread(filename_feb,'siv_S2026');

siv_A2010_jul = ncread(filename_jul,'siv_A2010');
siv_J2022_jul = ncread(filename_jul,'siv_J2021');
siv_Study_jul = ncread(filename_jul,'siv_S2026');

% Mean pan-Arctic sea-ice volume.
%
% Sea-ice volume is summed over the spatial grid separately for each year
% and then averaged over 2011-2023.
%
% The NetCDF sea-ice-volume fields are in m^3. Multiplication by 1e-12
% converts the result to 10^3 km^3.

V_A2010_feb = mean(sum(siv_A2010_feb,[1 2],'omitnan'),3,'omitnan') * 1e-12;
V_J2022_feb = mean(sum(siv_J2022_feb,[1 2],'omitnan'),3,'omitnan') * 1e-12;
V_Study_feb = mean(sum(siv_Study_feb,[1 2],'omitnan'),3,'omitnan') * 1e-12;

V_A2010_jul = mean(sum(siv_A2010_jul,[1 2],'omitnan'),3,'omitnan') * 1e-12;
V_J2022_jul = mean(sum(siv_J2022_jul,[1 2],'omitnan'),3,'omitnan') * 1e-12;
V_Study_jul = mean(sum(siv_Study_jul,[1 2],'omitnan'),3,'omitnan') * 1e-12;

% Relative sea-ice-volume differences.

dV_Study_A2010_feb = 100 * ...
    (V_Study_feb - V_A2010_feb) / V_A2010_feb;

dV_Study_A2010_jul = 100 * ...
    (V_Study_jul - V_A2010_jul) / V_A2010_jul;

dV_J2022_A2010_feb = 100 * ...
    (V_J2022_feb - V_A2010_feb) / V_A2010_feb;

dV_J2022_Study_feb = 100 * ...
    (V_J2022_feb - V_Study_feb) / V_Study_feb;

dV_J2022_A2010_jul = 100 * ...
    (V_J2022_jul - V_A2010_jul) / V_A2010_jul;

dV_J2022_Study_jul = 100 * ...
    (V_J2022_jul - V_Study_jul) / V_Study_jul;

% Climatological sea-ice thickness differences.
%
% Positive values mean that the first parameterization gives thicker ice.

dSIT_Study_A2010_feb = ...
    mean(sit_Study_feb - sit_A2010_feb,3,'omitnan');

dSIT_J2022_Study_feb = ...
    mean(sit_J2022_feb - sit_Study_feb,3,'omitnan');

dSIT_Study_A2010_jul = ...
    mean(sit_Study_jul - sit_A2010_jul,3,'omitnan');

dSIT_J2022_Study_jul = ...
    mean(sit_J2022_jul - sit_Study_jul,3,'omitnan');

% Apply the broad geographic and sea-ice-concentration mask used for the
% climatological thickness maps:
%   - exclude grid cells north of 88.5 degrees N,
%   - exclude cells with mean sea-ice concentration below 15%.
%
% Exact zero thickness differences are retained because zero is a valid
% difference between parameterizations.

mask_feb = grid_lat <= 88.5 & ...
    mean(sic_feb,3,'omitnan') >= 15;

mask_jul = grid_lat <= 88.5 & ...
    mean(sic_jul,3,'omitnan') >= 15;

x_Study_A2010_feb = ...
    dSIT_Study_A2010_feb(mask_feb & isfinite(dSIT_Study_A2010_feb));

x_J2022_Study_feb = ...
    dSIT_J2022_Study_feb(mask_feb & isfinite(dSIT_J2022_Study_feb));

x_Study_A2010_jul = ...
    dSIT_Study_A2010_jul(mask_jul & isfinite(dSIT_Study_A2010_jul));

x_J2022_Study_jul = ...
    dSIT_J2022_Study_jul(mask_jul & isfinite(dSIT_J2022_Study_jul));

% Percentiles are included because single-cell extrema can be sensitive to
% isolated grid cells. The 95th and 99th percentiles therefore provide
% additional context for the regional thickness-difference statements.

p = [5 50 95 99];

P_Study_A2010_feb = prctile(x_Study_A2010_feb,p);
P_J2022_Study_feb = prctile(x_J2022_Study_feb,p);

P_Study_A2010_jul = prctile(x_Study_A2010_jul,p);
P_J2022_Study_jul = prctile(x_J2022_Study_jul,p);

% Print values used to check the CryoSat-2 manuscript text.

fprintf('\n')
fprintf('============================================================\n')
fprintf('CRYOSAT-2 MANUSCRIPT VALUE CHECKS\n')
fprintf('============================================================\n')

fprintf('\nMean pan-Arctic sea-ice volume [10^3 km^3]\n')
fprintf('                         February       July\n')
fprintf('A2010                    %8.3f     %8.3f\n', ...
    V_A2010_feb,V_A2010_jul)
fprintf('J2022                    %8.3f     %8.3f\n', ...
    V_J2022_feb,V_J2022_jul)
fprintf('This Study               %8.3f     %8.3f\n', ...
    V_Study_feb,V_Study_jul)

fprintf('\nVolume differences [%%]\n')
fprintf('This Study vs A2010, February: %+6.2f %%\n', ...
    dV_Study_A2010_feb)
fprintf('This Study vs A2010, July:     %+6.2f %%\n', ...
    dV_Study_A2010_jul)
fprintf('J2022 vs A2010, February:      %+6.2f %%\n', ...
    dV_J2022_A2010_feb)
fprintf('J2022 vs This Study, February: %+6.2f %%\n', ...
    dV_J2022_Study_feb)
fprintf('J2022 vs A2010, July:          %+6.2f %%\n', ...
    dV_J2022_A2010_jul)
fprintf('J2022 vs This Study, July:     %+6.2f %%\n', ...
    dV_J2022_Study_jul)

fprintf('\nClimatological SIT difference ranges [cm]\n')
fprintf('This Study - A2010, February: %+.1f to %+.1f\n', ...
    100*min(x_Study_A2010_feb),100*max(x_Study_A2010_feb))
fprintf('J2022 - This Study, February: %+.1f to %+.1f\n', ...
    100*min(x_J2022_Study_feb),100*max(x_J2022_Study_feb))
fprintf('This Study - A2010, July:     %+.1f to %+.1f\n', ...
    100*min(x_Study_A2010_jul),100*max(x_Study_A2010_jul))
fprintf('J2022 - This Study, July:     %+.1f to %+.1f\n', ...
    100*min(x_J2022_Study_jul),100*max(x_J2022_Study_jul))

fprintf('\nSelected SIT-difference percentiles [cm]\n')
fprintf('                              5th    50th    95th    99th\n')
fprintf('This Study - A2010 Feb:    %6.1f %7.1f %7.1f %7.1f\n', ...
    100*P_Study_A2010_feb)
fprintf('J2022 - This Study Feb:    %6.1f %7.1f %7.1f %7.1f\n', ...
    100*P_J2022_Study_feb)
fprintf('This Study - A2010 Jul:    %6.1f %7.1f %7.1f %7.1f\n', ...
    100*P_Study_A2010_jul)
fprintf('J2022 - This Study Jul:    %6.1f %7.1f %7.1f %7.1f\n', ...
    100*P_J2022_Study_jul)

fprintf('============================================================\n')