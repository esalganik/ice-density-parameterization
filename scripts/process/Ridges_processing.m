function Ridges_processing(repoRoot)
%
% RIDGES_PROCESSING
%
% Imports ridge-core density and temperature measurements from MOSAiC
% ridge measurements, interpolates temperature profiles onto density-core
% section depths, computes in-situ density estimates from laboratory
% density, salinity, and temperature, and exports one summary record per
% ridge core.
%
% Input:
%   data/raw/ridges/**/*.xlsx
%
% Output:
%   data/processed/Summary_Ridges.mat
%

close all

rawDir = fullfile(repoRoot, 'data', 'raw', 'ridges');
processedDir = fullfile(repoRoot, 'data', 'processed');

if ~exist(processedDir, 'dir')
    mkdir(processedDir)
end

outputFile = fullfile(processedDir, 'Summary_Ridges.mat');

% Find all ridge-core workbooks and ignore temporary Excel files.
files = dir(fullfile(rawDir,'**','*.xlsx'));
files = files(~startsWith({files.name},'~$'));

cores = struct([]);

% Read ridge-core metadata, density sections, and temperature profiles from each workbook.
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
    if ~isnan(latDeg)
        if isnan(latMin)
            lat = latDeg;
        else
            lat = latDeg + latMin/60;
        end
    else
        lat = NaN;
    end

    if ~isnan(lonDeg)
        if isnan(lonMin)
            lon = lonDeg;
        else
            lon = lonDeg + lonMin/60;
        end
    else
        lon = NaN;
    end

    % Extract section-level density, salinity, and laboratory-temperature data.
    rowsDC = 23:size(dc,1);
    top = convcol(dc(rowsDC,1));
    bottom = convcol(dc(rowsDC,2));
    density = convcol(dc(rowsDC,8));
    salinity = convcol(dc(rowsDC,10));
    tempLab_section = -convcol(dc(rowsDC,12));

    keep = ~(isnan(top) & isnan(bottom) & isnan(density) & ...
        isnan(salinity) & isnan(tempLab_section));

    top = top(keep);
    bottom = bottom(keep);
    density = density(keep);
    salinity = salinity(keep);
    tempLab_section = tempLab_section(keep);

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
        temperature = min(temperature,-0.1);
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
    cores(k).temp_lab_section = tempLab_section;
    cores(k).snow_thickness = snowThickness;
    cores(k).ice_thickness = iceThickness;
    cores(k).ice_draft = iceDraft;
    cores(k).temperature_depth = depth;
    cores(k).temperature = temperature;
    cores(k).core_length = coreLength;
end

% Interpolate temperature profiles to density-section depths and compute
% in-situ density estimates for each ridge core.
for k = 1:numel(cores)
    zrho = (cores(k).density_top + cores(k).density_bottom) / 2;
    rho = cores(k).density;
    S = cores(k).salinity;
    zT = cores(k).temperature_depth;
    Tprof = cores(k).temperature;
    Tlab = cores(k).temp_lab_section;

    if isempty(Tlab) || all(isnan(Tlab))
        Tlab = cores(k).temp_lab * ones(size(rho));
    end

    T_rho = nan(size(rho));
    vb = nan(size(rho));
    vg = nan(size(rho));
    rho_si = nan(size(rho));

    if all(isnan(zT)) || all(isnan(Tprof))
        cores(k).temperature_rho = T_rho;
        cores(k).brine_volume = vb;
        cores(k).gas_volume = vg;
        cores(k).density_insitu = rho_si;
        cores(k).temperature_bulk = NaN;
        cores(k).density_insitu_bulk = NaN;
        continue
    end

    idxT = ~isnan(zT) & ~isnan(Tprof);
    idxR = ~isnan(zrho) & ~isnan(rho) & ~isnan(S) & ~isnan(Tlab);

    if sum(idxT) < 2 || sum(idxR) == 0
        cores(k).temperature_rho = T_rho;
        cores(k).brine_volume = vb;
        cores(k).gas_volume = vg;
        cores(k).density_insitu = rho_si;
        cores(k).temperature_bulk = NaN;
        cores(k).density_insitu_bulk = NaN;
        continue
    end

    zTok = zT(idxT);
    Tok = min(-0.1,Tprof(idxT));

    [zTuniq,~,ic] = unique(zTok);
    Tuniq = accumarray(ic,Tok,[],@mean);

    if numel(zTuniq) < 2
        cores(k).temperature_rho = T_rho;
        cores(k).brine_volume = vb;
        cores(k).gas_volume = vg;
        cores(k).density_insitu = rho_si;
        cores(k).temperature_bulk = NaN;
        cores(k).density_insitu_bulk = NaN;
        continue
    end

    % Interpolate in-situ temperatures from temperature-core measurements to the density-core section midpoints.
    T_rho(idxR) = interp1(zTuniq,Tuniq,zrho(idxR),'linear','extrap');

    TlabR = Tlab(idxR);
    Ti = T_rho(idxR);
    Si = S(idxR);
    rho_lab = rho(idxR);

    % Estimate brine volume, gas volume, and in-situ density from laboratory
    % density, salinity, laboratory temperature, and interpolated in-situ
    % temperature.
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

    vg_i = max(0,1 - (1 - vg_pr) .* A);

    rho_i = ...
        (1 - vg_i) .* rhoi_i .* F1_i ./ ...
        (F1_i - rhoi_i .* Si/1000 .* F2_i);

    rho_i(isnan(vb_calc)) = NaN;

    vb(idxR) = vb_export;
    vg(idxR) = vg_i;
    rho_si(idxR) = rho_i;

    cores(k).temperature_rho = T_rho;
    cores(k).brine_volume = vb;
    cores(k).gas_volume = vg;
    cores(k).density_insitu = rho_si;
    cores(k).temperature_bulk = mean(T_rho,'omitnan');
    cores(k).density_insitu_bulk = mean(rho_si,'omitnan');
end

n = numel(cores);

rho_ridge = nan(n,1);
rho_si_ridge = nan(n,1);
hi_ridge = nan(n,1);
hs_ridge = nan(n,1);
Ti_ridge = nan(n,1);
t_ridge = NaT(n,1);
lat_ridge = nan(n,1);
lon_ridge = nan(n,1);
dataset_ridge = strings(n,1);
file_ridge = strings(n,1);

for k = 1:n
    file_ridge(k) = string(cores(k).file);
    dataset_ridge(k) = string(cores(k).source);

    % Average section-level properties to one representative value per ridge core.
    rho_ridge(k) = mean(cores(k).density,'omitnan');
    rho_si_ridge(k) = mean(cores(k).density_insitu,'omitnan');

    hi_ridge(k) = cores(k).ice_thickness / 100;
    hs_ridge(k) = cores(k).snow_thickness / 100;
    Ti_ridge(k) = mean(cores(k).temperature_rho,'omitnan');

    t_ridge(k) = cores(k).date;
    lat_ridge(k) = cores(k).lat;
    lon_ridge(k) = cores(k).lon;
end

% Convert processed ridge-core means to the common repository schema.
Ridge_core_summary = table( ...
    file_ridge, ...
    dataset_ridge, ...
    t_ridge, ...
    rho_ridge, ...
    rho_si_ridge, ...
    hi_ridge, ...
    hs_ridge, ...
    Ti_ridge, ...
    lat_ridge, ...
    lon_ridge, ...
'VariableNames', { ...
    'file', ...
    'dataset', ...
    'date', ...
    'density_lab_kgm3', ...
    'density_insitu_kgm3', ...
    'ice_thickness_m', ...
    'snow_thickness_m', ...
    'temperature_C', ...
    'lat_degN', ...
    'lon_degE'});

Ridge_core_summary = sortrows(Ridge_core_summary,{'dataset','date'});

if ismember('file',Ridge_core_summary.Properties.VariableNames) && ...
        ~ismember('core',Ridge_core_summary.Properties.VariableNames)
    Ridge_core_summary.core = Ridge_core_summary.file;
    Ridge_core_summary = movevars(Ridge_core_summary,'core','After','dataset');
end

Summary_Ridges = Ridge_core_summary;

fprintf('Imported ridge dataset.\n')
fprintf('Processed %d ice cores.\n',height(Summary_Ridges))
fprintf('Output saved to:\n%s\n',outputFile)

save(outputFile,'Summary_Ridges')

%% Helper functions

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