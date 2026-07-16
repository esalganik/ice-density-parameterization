function MOSAiC_SYI_processing(repoRoot)
%
% MOSAIC_SYI_PROCESSING
%
% Imports MOSAiC second-year ice core density data, computes in-situ density
% estimates from laboratory density, salinity, and in-situ temperature, and
% averages section-level measurements to one record per core/date.
%
% Input:
%   data/raw/MOSAiC_MCS_SYI_sea_ice.tab
%
% Output:
%   data/processed/Summary_MOSAIC_SYI.mat
%

close all

rawDir = fullfile(repoRoot, 'data', 'raw');
processedDir = fullfile(repoRoot, 'data', 'processed');

if ~exist(processedDir, 'dir')
    mkdir(processedDir)
end

inputFile = fullfile(rawDir, ...
    'MOSAiC_MCS_SYI_sea_ice.tab');

outputFile = fullfile(processedDir, ...
    'Summary_MOSAIC_SYI.mat');

% Remove PANGAEA metadata header and read the tabular data section only.
txt = fileread(inputFile);
endHeader = strfind(txt,'*/');
dataTxt = strtrim(txt(endHeader(1)+2:end));

tmpFile = [tempname '.tab'];
fid = fopen(tmpFile,'w');
fwrite(fid,dataTxt);
fclose(fid);

opts = detectImportOptions(tmpFile, ...
    'FileType','text', ...
    'Delimiter','\t', ...
    'VariableNamingRule','preserve');

Tfull = readtable(tmpFile,opts);
delete(tmpFile)

% Select variables needed for density processing and repository summary output.
cols = [2 3 4 5 8 9 14 19 20 22 23];
T = Tfull(:,cols);

T.Properties.VariableNames = { ...
    'DateTime', ...
    'Latitude', ...
    'Longitude', ...
    'SeaIceThickness_m', ...
    'DepthTop_m', ...
    'DepthBottom_m', ...
    'SeaIceSalinity', ...
    'Temperature_C', ...
    'Density_lab_kgm3', ...
    'SnowThickness_m', ...
    'Comment'};

T.DateTime = datetime(T.DateTime, ...
    'InputFormat','yyyy-MM-dd''T''HH:mm');

% Laboratory density measurements are converted using a fixed reference laboratory temperature because no sample-specific 
% laboratory temperature is provided in this dataset.
Tlab_fixed = -15;

S = T.SeaIceSalinity;
rho = T.Density_lab_kgm3;
% Limit in-situ ice temperature to below the freezing point to avoid invalid brine-volume calculations near 0 deg C.
Tinsitu = min(-0.1,T.Temperature_C);

T.Temperature_C_limited = Tinsitu;
T.Tlab_C = Tlab_fixed * ones(height(T),1);

T.BrineVolume = nan(height(T),1);
T.GasVolume_lab = nan(height(T),1);
T.GasVolume_connected = nan(height(T),1);
T.GasVolume_disconnected = nan(height(T),1);
T.Density_insitu_connected_kgm3 = nan(height(T),1);
T.Density_insitu_disconnected_kgm3 = nan(height(T),1);

idx = ~isnan(S) & ~isnan(rho) & ~isnan(Tinsitu);

Si = S(idx);
rho_lab = rho(idx);
Tlab_i = Tlab_fixed * ones(sum(idx),1);
Ti_i = Tinsitu(idx);

% Estimate brine volume, gas volume, and in-situ density from laboratory
% density, salinity, fixed laboratory temperature, and in-situ temperature.
[F1_pr,F2_pr] = F1F2_seaice(Tlab_i);
[F1_i,F2_i] = F1F2_seaice(Ti_i);

rhoi_pr = 917 - 0.1403*Tlab_i;
rhoi_i  = 917 - 0.1403*Ti_i;

vb_pr = rho_lab .* Si ./ F1_pr;

vg_pr = max(0, ...
    1 - rho_lab .* ...
    (F1_pr - rhoi_pr .* Si/1000 .* F2_pr) ./ ...
    (rhoi_pr .* F1_pr));

F3_pr = rhoi_pr .* Si/1000 ./ ...
    (F1_pr - rhoi_pr .* Si/1000 .* F2_pr);

F3_i = rhoi_i .* Si/1000 ./ ...
    (F1_i - rhoi_i .* Si/1000 .* F2_i);

vb_raw = vb_pr .* F1_pr ./ F1_i / 1000;

vb_export = vb_raw;
vb_export(vb_export > 0.4 | vb_export < 0) = NaN;

vb_calc = vb_raw;
vb_calc(vb_calc > 0.4 | vb_calc < 0) = 0.4;

% Connected-pore case: allow gas volume to change when laboratory density
% is converted to in-situ temperature conditions.
A = ...
    (rhoi_i ./ rhoi_pr) .* ...
    (F3_pr ./ F3_i) .* ...
    (F1_pr ./ F1_i);

idx0 = Si == 0;
A(idx0) = 1;

vg_connected = max(0, ...
    1 - (1 - vg_pr) .* A);

% Disconnected-pore case: preserve laboratory gas volume during conversion
% to in-situ temperature conditions.
vg_disconnected = vg_pr;

rho_si_connected = ...
    (1 - vg_connected) .* rhoi_i .* F1_i ./ ...
    (F1_i - rhoi_i .* Si/1000 .* F2_i);

rho_si_disconnected = ...
    (1 - vg_disconnected) .* rhoi_i .* F1_i ./ ...
    (F1_i - rhoi_i .* Si/1000 .* F2_i);

rho_si_connected(isnan(vb_calc)) = NaN;
rho_si_disconnected(isnan(vb_calc)) = NaN;

T.BrineVolume(idx) = vb_export;
T.GasVolume_lab(idx) = vg_pr;
T.GasVolume_connected(idx) = vg_connected;
T.GasVolume_disconnected(idx) = vg_disconnected;
T.Density_insitu_connected_kgm3(idx) = rho_si_connected;
T.Density_insitu_disconnected_kgm3(idx) = rho_si_disconnected;

% Exclude false-bottom sections and missing laboratory densities before
% averaging section-level measurements to core/date means.
idxFalseBottom = contains(lower(string(T.Comment)),'false bottom');
idxNaNDensity = isnan(T.Density_lab_kgm3);

TavgInput = T(~idxFalseBottom & ~idxNaNDensity,:);

numericVars = varfun(@isnumeric,TavgInput,'OutputFormat','uniform');

excludeFromAverages = ismember(TavgInput.Properties.VariableNames, ...
    {'DepthTop_m','DepthBottom_m'});

numNames = TavgInput.Properties.VariableNames(numericVars & ~excludeFromAverages);

% Average numeric section-level variables to one representative record per core/date.
Tavg = groupsummary(TavgInput,'DateTime','mean',numNames);

if ismember('GroupCount',Tavg.Properties.VariableNames)
    Tavg.GroupCount = [];
end

Tavg.Properties.VariableNames = erase(Tavg.Properties.VariableNames,'mean_');

fprintf('Imported MOSAiC SYI dataset.\n')
fprintf('Processed %d ice cores.\n', ...
    height(Tavg))
fprintf('Output saved to:\n%s\n', ...
    outputFile)

dataset = repmat("MOSAiC SYI",height(Tavg),1);
core = string(Tavg.DateTime);
date = Tavg.DateTime;
ice_age = repmat("SYI",height(Tavg),1);

temperature_C = Tavg.Temperature_C_limited;
density_lab_kgm3 = Tavg.Density_lab_kgm3;
density_insitu_connected_kgm3 = Tavg.Density_insitu_connected_kgm3;
density_insitu_disconnected_kgm3 = Tavg.Density_insitu_disconnected_kgm3;
ice_thickness_m = Tavg.SeaIceThickness_m;
snow_thickness_m = Tavg.SnowThickness_m;
lat_degN = Tavg.Latitude;
lon_degE = Tavg.Longitude;

Summary_MOSAIC_SYI = table( ...
    dataset, ...
    core, ...
    date, ...
    ice_age, ...
    temperature_C, ...
    density_lab_kgm3, ...
    density_insitu_connected_kgm3, ...
    density_insitu_disconnected_kgm3, ...
    ice_thickness_m, ...
    snow_thickness_m, ...
    lat_degN, ...
    lon_degE, ...
    'VariableNames',{ ...
    'dataset', ...
    'core', ...
    'date', ...
    'ice_age', ...
    'temperature_C', ...
    'density_lab_kgm3', ...
    'density_insitu_connected_kgm3', ...
    'density_insitu_disconnected_kgm3', ...
    'ice_thickness_m', ...
    'snow_thickness_m', ...
    'lat_degN', ...
    'lon_degE'});

save(outputFile,'Summary_MOSAIC_SYI')

%% Helpers
% Compute empirical sea-ice brine-volume coefficients F1 and F2 for the
% relevant temperature regime.
function [F1,F2] = F1F2_seaice(T)

F1 = -4.732 - 22.45*T - 0.6397*T.^2 - 0.01074*T.^3;
F2 = 8.903e-2 - 1.763e-2*T - 5.33e-4*T.^2 - 8.801e-6*T.^3;

idxWarm = T > -2;

F1(idxWarm) = -4.1221e-2 - 18.407*T(idxWarm) + ...
    0.58402*T(idxWarm).^2 + 0.21454*T(idxWarm).^3;

F2(idxWarm) = 9.0312e-2 - 0.016111*T(idxWarm) + ...
    1.2291e-4*T(idxWarm).^2 + 1.3603e-4*T(idxWarm).^3;

idxCold = T < -22.9;

F1(idxCold) = 9899 + 1309*T(idxCold) + ...
    55.27*T(idxCold).^2 + 0.7160*T(idxCold).^3;

F2(idxCold) = 8.547 + 1.089*T(idxCold) + ...
    0.04518*T(idxCold).^2 + 5.819e-4*T(idxCold).^3;

end
end