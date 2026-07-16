function Svalbard_processing(repoRoot)
%
% SVALBARD_PROCESSING
%
% Imports sea-ice density and temperature measurements from Svalbard and
% N-ICE2015 ice-core campaigns, interpolates temperature profiles onto
% density-core section depths, computes in-situ density estimates from
% laboratory density, salinity, and temperature, and exports one summary
% record per ice core.
%
% Input:
%   data/raw/Svalbard/**/*.xlsx
%
% Output:
%   data/processed/Summary_Svalbard.mat
%

close all

rawDir = fullfile(repoRoot, 'data', 'raw', 'Svalbard');
processedDir = fullfile(repoRoot, 'data', 'processed');

if ~exist(processedDir, 'dir')
    mkdir(processedDir)
end

outputFile = fullfile(processedDir, ...
    'Summary_Svalbard.mat');

rootFolder = rawDir;

folder1 = fullfile(rootFolder,'Kongsfjorden_density_cores_for_Evgenii');
folder2 = fullfile(rootFolder,'Density_cores_KF2018_NICE');

% Collect all Svalbard and N-ICE2015 workbooks.
files = [ ...
    dir(fullfile(folder1,'*.xlsx')); ...
    dir(fullfile(folder2,'*.xlsx'))];

cores = struct([]);

% Read metadata, density-core measurements, and temperature profiles from each workbook.
for k = 1:numel(files)
    file = fullfile(files(k).folder,files(k).name);

    dc = readcell(file,'Sheet','Density Core');
    [~,sheets] = xlsfinfo(file);
    hasTemp = any(strcmp(sheets,'Temperature Core'));

    so = readcell(file,'Sheet','Station Overview');

    latDeg = convcell(so{20,3});
    latMin = convcell(so{20,4});
    lonDeg = convcell(so{21,3});
    lonMin = convcell(so{21,4});

    % Convert station coordinates from degrees and minutes to decimal degrees.
    lat = NaN;
    lon = NaN;

    if ~isnan(latDeg) && ~isnan(latMin)
        lat = latDeg + latMin/60;
    end

    if ~isnan(lonDeg) && ~isnan(lonMin)
        lon = lonDeg + lonMin/60;
    end

    % Extract section-level density and salinity measurements.
    rowsDC = 23:size(dc,1);
    top = convcol(dc(rowsDC,1));
    bottom = convcol(dc(rowsDC,2));
    density = 1000 * convcol(dc(rowsDC,8));
    salinity = convcol(dc(rowsDC,10));

    keep = ~(isnan(top) & isnan(bottom) & isnan(density) & isnan(salinity));

    top = top(keep);
    bottom = bottom(keep);
    density = density(keep);
    salinity = salinity(keep);

    tempLab = convcell(dc{7,7});
    snowThickness = convcell(dc{10,2});
    iceThickness = convcell(dc{11,2});
    iceDraft = convcell(dc{12,2});

    dateCore = parsedate(dc{4,2});

    % Extract temperature-profile measurements when available.
    if hasTemp
        tc = readcell(file,'Sheet','Temperature Core');
        coreLength = convcell(tc{18,5});
        rowsTC = 23:size(tc,1);

        depthRaw = tc(rowsTC,1);
        temperature = convcol(tc(rowsTC,2));
        depth = nan(size(temperature));

        for i = 1:numel(depthRaw)
            v = depthRaw{i};

            if isempty(v)
                depth(i) = NaN;
            elseif isnumeric(v)
                depth(i) = double(v);
            else
                s = lower(strtrim(string(v)));

                if s == "top"
                    depth(i) = 0;
                elseif s == "bottom"
                    depth(i) = coreLength;
                else
                    depth(i) = str2double(regexprep(s,'[^0-9\.\-]',''));
                end
            end
        end

        keepT = ~(isnan(depth) & isnan(temperature));
        depth = depth(keepT);
        temperature = temperature(keepT);
    else
        depth = NaN;
        temperature = NaN;
        coreLength = NaN;
    end

    [~,sourceName] = fileparts(files(k).folder);

    cores(k).file = files(k).name;
    cores(k).source = sourceName;
    cores(k).date = dateCore;
    cores(k).lat = lat;
    cores(k).lon = lon;
    cores(k).density_top = top;
    cores(k).density_bottom = bottom;
    cores(k).density = density;
    cores(k).salinity = salinity;
    cores(k).temp_lab = tempLab;
    cores(k).snow_thickness = snowThickness;
    cores(k).ice_thickness = iceThickness;
    cores(k).ice_draft = iceDraft;
    cores(k).temperature_depth = depth;
    cores(k).temperature = temperature;
    cores(k).core_length = coreLength;
end

% Interpolate temperature profiles to density-section depths and compute in-situ density estimates for each ice core.
for k = 1:numel(cores)

    zrho = (cores(k).density_top + cores(k).density_bottom) / 2;
    rho = cores(k).density;
    S = cores(k).salinity;
    zT = cores(k).temperature_depth;
    Tprof = cores(k).temperature;
    Tlab = cores(k).temp_lab * ones(size(rho));

    T_rho = nan(size(rho));
    vb = nan(size(rho));
    vg_lab = nan(size(rho));
    vg_connected = nan(size(rho));
    vg_disconnected = nan(size(rho));
    rho_si_connected = nan(size(rho));
    rho_si_disconnected = nan(size(rho));

    idxT = ~isnan(zT) & ~isnan(Tprof);
    idxR = ~isnan(zrho) & ~isnan(rho) & ~isnan(S) & ~isnan(Tlab);

    if sum(idxT) < 2 || sum(idxR) == 0
        cores(k).temperature_rho = T_rho;
        cores(k).brine_volume = vb;
        cores(k).gas_volume_lab = vg_lab;
        cores(k).gas_volume_connected = vg_connected;
        cores(k).gas_volume_disconnected = vg_disconnected;
        cores(k).density_insitu_connected = rho_si_connected;
        cores(k).density_insitu_disconnected = rho_si_disconnected;
        continue
    end

    zTok = zT(idxT);
    Tok = min(-0.1,Tprof(idxT));

    [zTuniq,~,ic] = unique(zTok);
    Tuniq = accumarray(ic,Tok,[],@mean);

    if numel(zTuniq) < 2
        cores(k).temperature_rho = T_rho;
        cores(k).brine_volume = vb;
        cores(k).gas_volume_lab = vg_lab;
        cores(k).gas_volume_connected = vg_connected;
        cores(k).gas_volume_disconnected = vg_disconnected;
        cores(k).density_insitu_connected = rho_si_connected;
        cores(k).density_insitu_disconnected = rho_si_disconnected;
        continue
    end

    % Interpolate in-situ temperatures from temperature-core measurements to density-core section midpoints.
    T_rho(idxR) = interp1(zTuniq,Tuniq,zrho(idxR),'linear','extrap');

    TlabR = Tlab(idxR);
    Ti = T_rho(idxR);
    Si = S(idxR);
    rho_lab = rho(idxR);

    % Estimate brine volume, gas volume, and in-situ density from laboratory
    % density, salinity, laboratory temperature, and interpolated in-situ temperature.

    [F1_pr,F2_pr] = F1F2_seaice(TlabR);
    [F1_i,F2_i] = F1F2_seaice(Ti);

    rhoi_pr = 917 - 0.1403*TlabR;
    rhoi_i = 917 - 0.1403*Ti;

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

    % Connected-pore case: allow gas volume to change when laboratory density is converted to in-situ temperature conditions.
    A = ...
        (rhoi_i ./ rhoi_pr) .* ...
        (F3_pr ./ F3_i) .* ...
        (F1_pr ./ F1_i);

    idx0 = Si == 0;
    A(idx0) = 1;

    vg_conn = max(0, ...
        1 - (1 - vg_pr) .* A);

    % Disconnected-pore case: preserve laboratory gas volume during conversion to in-situ temperature conditions.
    vg_disc = vg_pr;

    rho_conn = ...
        (1 - vg_conn) .* rhoi_i .* F1_i ./ ...
        (F1_i - rhoi_i .* Si/1000 .* F2_i);

    rho_disc = ...
        (1 - vg_disc) .* rhoi_i .* F1_i ./ ...
        (F1_i - rhoi_i .* Si/1000 .* F2_i);

    rho_conn(isnan(vb_calc)) = NaN;
    rho_disc(isnan(vb_calc)) = NaN;

    vb(idxR) = vb_export;
    vg_lab(idxR) = vg_pr;
    vg_connected(idxR) = vg_conn;
    vg_disconnected(idxR) = vg_disc;
    rho_si_connected(idxR) = rho_conn;
    rho_si_disconnected(idxR) = rho_disc;

    cores(k).temperature_rho = T_rho;
    cores(k).brine_volume = vb;
    cores(k).gas_volume_lab = vg_lab;
    cores(k).gas_volume_connected = vg_connected;
    cores(k).gas_volume_disconnected = vg_disconnected;
    cores(k).density_insitu_connected = rho_si_connected;
    cores(k).density_insitu_disconnected = rho_si_disconnected;
    cores(k).temperature_bulk = mean(T_rho,'omitnan');
    cores(k).density_insitu_connected_bulk = mean(rho_si_connected,'omitnan');
    cores(k).density_insitu_disconnected_bulk = mean(rho_si_disconnected,'omitnan');

end

nCores = numel(cores);

file_sva = strings(nCores,1);
dataset_sva = strings(nCores,1);
coreDate = NaT(nCores,1);

avgDensity = nan(nCores,1);
avgDensity_insitu_connected = nan(nCores,1);
avgDensity_insitu_disconnected = nan(nCores,1);
iceThickness = nan(nCores,1);
snowDepth = nan(nCores,1);
avgTemp = nan(nCores,1);
lat_sva = nan(nCores,1);
lon_sva = nan(nCores,1);

% Average section-level properties to one representative value per ice core.
for k = 1:nCores
    file_sva(k) = string(cores(k).file);

    % Separate N-ICE2015 cores from the remaining Svalbard observations.
    if strcmp(cores(k).file,'N-ICE_Ice_Station_2015_02_12_with_density_100426.xlsx')
        dataset_sva(k) = "N-ICE2015";
    else
        dataset_sva(k) = "Svalbard";
    end

    coreDate(k) = cores(k).date;

    avgDensity(k) = mean(cores(k).density,'omitnan');
    avgDensity_insitu_connected(k) = mean(cores(k).density_insitu_connected,'omitnan');
    avgDensity_insitu_disconnected(k) = mean(cores(k).density_insitu_disconnected,'omitnan');

    iceThickness(k) = cores(k).ice_thickness;
    snowDepth(k) = cores(k).snow_thickness;
    avgTemp(k) = mean(cores(k).temperature_rho,'omitnan');

    lat_sva(k) = cores(k).lat;
    lon_sva(k) = cores(k).lon;
end

hi_sva = iceThickness / 100;
hs_sva = snowDepth / 100;
Ti_sva = avgTemp;
rho_sva = avgDensity;
ice_age = repmat("FYI",nCores,1);

% Convert processed core means to the common repository schema.
Svalbard_core_summary = table( ...
    file_sva, ...
    dataset_sva, ...
    coreDate, ...
    ice_age, ...
    rho_sva, ...
    avgDensity_insitu_connected, ...
    avgDensity_insitu_disconnected, ...
    hi_sva, ...
    hs_sva, ...
    Ti_sva, ...
    lat_sva, ...
    lon_sva, ...
'VariableNames', { ...
    'file', ...
    'dataset', ...
    'date', ...
    'ice_age', ...
    'density_lab_kgm3', ...
    'density_insitu_connected_kgm3', ...
    'density_insitu_disconnected_kgm3', ...
    'ice_thickness_m', ...
    'snow_depth_m', ...
    'temperature_C', ...
    'lat_degN', ...
    'lon_degE'});

Svalbard_core_summary = sortrows(Svalbard_core_summary,{'dataset','date'});

fprintf('Imported Svalbard dataset.\n')
fprintf('Processed %d ice cores.\n', ...
    height(Svalbard_core_summary))
fprintf('Output saved to:\n%s\n', ...
    outputFile)

save(outputFile,'Svalbard_core_summary')

%% Helpers

% F1F2_SEAICE
%
% Computes the empirical sea-ice coefficients F1 and F2 used in the Cox
% and Weeks brine-volume formulation as a function of ice temperature.
%
% Input:
%   T  - Sea-ice temperature (°C)
%
% Output:
%   F1 - Temperature-dependent coefficient
%   F2 - Temperature-dependent coefficient
%
% Notes:
%   Piecewise polynomial fits are used following Cox and Weeks (1983):
%       -22.9 °C <= T <= -2 °C
%       T > -2 °C
%       T < -22.9 °C
%
%   These coefficients are used to estimate brine volume, gas volume, and
%   temperature-corrected in-situ sea-ice density from laboratory density
%   measurements.
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

function x = convcol(c)
x = nan(size(c));
for ii = 1:numel(c)
    x(ii) = convcell(c{ii});
end
end

function x = convcell(v)
if isempty(v)
    x = NaN;
elseif isnumeric(v)
    x = double(v);
elseif islogical(v)
    x = double(v);
elseif ischar(v) || isstring(v)
    x = str2double(string(v));
else
    x = NaN;
end
end

function dateCore = parsedate(dateRaw)
if isdatetime(dateRaw)
    dateCore = dateRaw;
elseif isnumeric(dateRaw)
    dateCore = datetime(dateRaw,'ConvertFrom','excel');
elseif ischar(dateRaw) || isstring(dateRaw)
    s = strtrim(string(dateRaw));

    try
        dateCore = datetime(s,'InputFormat','dd-MM-yy');
    catch
        try
            dateCore = datetime(s,'InputFormat','dd.MM.yyyy');
        catch
            try
                dateCore = datetime(s);
            catch
                dateCore = NaT;
            end
        end
    end
else
    dateCore = NaT;
end
end

end