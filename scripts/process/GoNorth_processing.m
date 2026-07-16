function GoNorth_processing(repoRoot)
%
% GONORTH_PROCESSING
%
% Imports GoNorth sea-ice core density data, interpolates temperature
% profiles onto density-core section depths, computes in-situ density
% estimates from laboratory density, salinity, and temperature, and averages
% section-level measurements to one record per ice core.
%
% Input:
%   data/raw/GoNorth1_icecores.tab
%
% Output:
%   data/processed/Summary_GoNorth.mat
%

close all

rawDir = fullfile(repoRoot, 'data', 'raw');
processedDir = fullfile(repoRoot, 'data', 'processed');

if ~exist(processedDir, 'dir')
    mkdir(processedDir)
end

inputFile = fullfile(rawDir, ...
    'GoNorth1_icecores.tab');

outputFile = fullfile(processedDir, ...
    'Summary_GoNorth.mat');

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

Tfull.Properties.VariableNames = { ...
    'Event','Station','DateTime','Latitude','Longitude', ...
    'SeaIceThickness_m','SeaIceDraft_m','SnowThickness_m', ...
    'DistX_m','DistY_m','IceAge','Core', ...
    'DepthTop_m','DepthBottom_m','SeaIceSalinity', ...
    'Temperature_C','Tlab_C','Density_lab_kgm3', ...
    'BrineVolume_file','GasVolume_file'};

Tfull.DateTime = datetime(Tfull.DateTime, ...
    'InputFormat','yyyy-MM-dd''T''HH:mm');

Tfull.IceAge = upper(string(Tfull.IceAge));
Tfull.Core = string(Tfull.Core);
Tfull.Station = string(Tfull.Station);

% Separate density cores from temperature profiles. Temperature profiles are later interpolated onto the density-core section depths.
Tden = Tfull(startsWith(Tfull.Core,"DEN"),:);
Ttemp = Tfull(Tfull.Core == "TEMP",:);

Tden.Temperature_C_interp = nan(height(Tden),1);
Tden.BrineVolume = nan(height(Tden),1);
Tden.GasVolume_lab = nan(height(Tden),1);
Tden.GasVolume_connected = nan(height(Tden),1);
Tden.GasVolume_disconnected = nan(height(Tden),1);
Tden.Density_insitu_connected_kgm3 = nan(height(Tden),1);
Tden.Density_insitu_disconnected_kgm3 = nan(height(Tden),1);

% Interpolate the station temperature profile to each density-core section midpoint. Profile depths are scaled to the density-core thickness before interpolation.
coreKeys = unique(strcat(Tden.Station,"_",Tden.Core),'stable');

for k = 1:numel(coreKeys)
    idxCore = strcat(Tden.Station,"_",Tden.Core) == coreKeys(k);
    stationNow = Tden.Station(find(idxCore,1,'first'));

    idxTemp = Ttemp.Station == stationNow;

    zrho = mean([Tden.DepthTop_m(idxCore),Tden.DepthBottom_m(idxCore)],2);

    zT = Ttemp.DepthTop_m(idxTemp);
    temp = Ttemp.Temperature_C(idxTemp);
    temp = min(-0.1,temp);

    okT = ~isnan(zT) & ~isnan(temp);

    if sum(okT) >= 2
        zTok = zT(okT);
        tempok = temp(okT);

        zTok = zTok * max(zrho,[],'omitnan') / max(zTok,[],'omitnan');

        [zTuniq,~,ic] = unique(zTok);
        tempuniq = accumarray(ic,tempok,[],@mean);

        Tden.Temperature_C_interp(idxCore) = interp1( ...
            zTuniq,tempuniq,zrho,'linear','extrap');
    end
end

S = Tden.SeaIceSalinity;
rho = Tden.Density_lab_kgm3;
Tlab = Tden.Tlab_C;
Ti = Tden.Temperature_C_interp;

idx = ~isnan(S) & ~isnan(rho) & ~isnan(Tlab) & ~isnan(Ti);

Si = S(idx);
rho_lab = rho(idx);
Tlab_i = Tlab(idx);
Ti_i = Ti(idx);

% Estimate brine volume, gas volume, and in-situ density from laboratory
% density, salinity, laboratory temperature, and interpolated in-situ
% temperature.
[F1_pr,F2_pr] = F1F2_seaice(Tlab_i);
[F1_i,F2_i] = F1F2_seaice(Ti_i);

rhoi_pr = 917 - 0.1403*Tlab_i;
rhoi_i = 917 - 0.1403*Ti_i;

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
vb_export(vb_export > 0.6 | vb_export < 0) = NaN;

vb_calc = vb_raw;
vb_calc(vb_calc > 0.6 | vb_calc < 0) = 0.6;

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

% Disconnected-pore case: preserve laboratory gas volume during conversion to in-situ temperature conditions.
vg_disconnected = vg_pr;

rho_connected = ...
    (1 - vg_connected) .* rhoi_i .* F1_i ./ ...
    (F1_i - rhoi_i .* Si/1000 .* F2_i);

rho_disconnected = ...
    (1 - vg_disconnected) .* rhoi_i .* F1_i ./ ...
    (F1_i - rhoi_i .* Si/1000 .* F2_i);

rho_connected(isnan(vb_calc)) = NaN;
rho_disconnected(isnan(vb_calc)) = NaN;

Tden.BrineVolume(idx) = vb_export;
Tden.GasVolume_lab(idx) = vg_pr;
Tden.GasVolume_connected(idx) = vg_connected;
Tden.GasVolume_disconnected(idx) = vg_disconnected;
Tden.Density_insitu_connected_kgm3(idx) = rho_connected;
Tden.Density_insitu_disconnected_kgm3(idx) = rho_disconnected;

% Average density-section measurements to one representative record per station/core.
Tden.CoreID = strcat(Tden.Station,"_",Tden.Core);

numericVars = varfun(@isnumeric,Tden,'OutputFormat','uniform');

excludeVars = ismember(Tden.Properties.VariableNames, ...
    {'DepthTop_m','DepthBottom_m','SeaIceDraft_m', ...
    'DistX_m','DistY_m','BrineVolume_file','GasVolume_file'});

numNames = Tden.Properties.VariableNames(numericVars & ~excludeVars);

Tavg = groupsummary(Tden,'CoreID','mean',numNames);

if ismember('GroupCount',Tavg.Properties.VariableNames)
    Tavg.GroupCount = [];
end

Tavg.Properties.VariableNames = erase(Tavg.Properties.VariableNames,'mean_');

coreInfo = groupsummary(Tden(:,{'CoreID','DateTime'}),'CoreID','min','DateTime');
coreInfo.GroupCount = [];
coreInfo.Properties.VariableNames{'min_DateTime'} = 'DateTime';

iceAgeInfo = groupsummary(Tden(:,{'CoreID','IceAge'}),'CoreID', ...
    @(x) string(x(find(string(x) ~= "",1,'first'))),'IceAge');
iceAgeInfo.GroupCount = [];
iceAgeInfo.Properties.VariableNames{'fun1_IceAge'} = 'IceAge';

Tavg = outerjoin(coreInfo,Tavg,'Keys','CoreID','MergeKeys',true,'Type','right');
Tavg = outerjoin(iceAgeInfo,Tavg,'Keys','CoreID','MergeKeys',true,'Type','right');

dataset = repmat("GoNorth SYI",height(Tavg),1);
idxFYI = upper(string(Tavg.IceAge)) == "FYI";
dataset(idxFYI) = "GoNorth FYI";
core = Tavg.CoreID;
date = Tavg.DateTime;
ice_age = upper(string(Tavg.IceAge));
temperature_C = Tavg.Temperature_C_interp;
density_lab_kgm3 = Tavg.Density_lab_kgm3;
density_insitu_connected_kgm3 = Tavg.Density_insitu_connected_kgm3;
density_insitu_disconnected_kgm3 = Tavg.Density_insitu_disconnected_kgm3;
ice_thickness_m = Tavg.SeaIceThickness_m;
snow_thickness_m = Tavg.SnowThickness_m;
lat_degN = Tavg.Latitude;
lon_degE = Tavg.Longitude;

Summary_GoNorth = table( ...
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

fprintf('Imported GoNorth dataset.\n')
fprintf('Processed %d ice cores.\n', ...
    height(Summary_GoNorth))
fprintf('Output saved to:\n%s\n', ...
    outputFile)

save(outputFile,'Summary_GoNorth')

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