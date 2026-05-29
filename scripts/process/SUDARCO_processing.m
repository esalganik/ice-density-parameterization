function SUDARCO_processing(repoRoot)

close all

rawDir = fullfile(repoRoot, 'data', 'raw', 'SUDARCO');
processedDir = fullfile(repoRoot, 'data', 'processed');

if ~exist(processedDir, 'dir')
    mkdir(processedDir)
end

outputFile = fullfile(processedDir, ...
    'Summary_SUDARCO.mat');

dataFolder = rawDir;

files = dir(fullfile(dataFolder,'**','*.xlsx'));
files = files(~startsWith({files.name},'~$'));

cores = struct([]);

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

    lat = NaN;
    lon = NaN;

    if ~isnan(latDeg) && ~isnan(latMin)
        lat = latDeg + latMin/60;
    end

    if ~isnan(lonDeg) && ~isnan(lonMin)
        lon = lonDeg + lonMin/60;
    end

    rowsDC = 23:size(dc,1);
    top = convcol(dc(rowsDC,1));
    bottom = convcol(dc(rowsDC,2));
    density = 1000 * convcol(dc(rowsDC,8));
    salinity = convcol(dc(rowsDC,10));
    tempLab_section = -convcol(dc(rowsDC,12));

    keep = ~(isnan(top) & isnan(bottom) & isnan(density) & isnan(salinity) & isnan(tempLab_section));

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
    cores(k).temp_lab_section = tempLab_section;
    cores(k).snow_thickness = snowThickness;
    cores(k).ice_thickness = iceThickness;
    cores(k).ice_draft = iceDraft;
    cores(k).temperature_depth = depth;
    cores(k).temperature = temperature;
    cores(k).core_length = coreLength;
end

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
    vg_pr_out = nan(size(rho));
    vg_connected = nan(size(rho));
    vg_disconnected = nan(size(rho));
    rho_si_connected = nan(size(rho));
    rho_si_disconnected = nan(size(rho));

    if all(isnan(zT)) || all(isnan(Tprof))
        cores(k).temperature_rho = T_rho;
        cores(k).brine_volume = vb;
        cores(k).gas_volume_lab = vg_pr_out;
        cores(k).gas_volume_connected = vg_connected;
        cores(k).gas_volume_disconnected = vg_disconnected;
        cores(k).density_insitu_connected = rho_si_connected;
        cores(k).density_insitu_disconnected = rho_si_disconnected;
        cores(k).temperature_bulk = NaN;
        cores(k).density_insitu_connected_bulk = NaN;
        cores(k).density_insitu_disconnected_bulk = NaN;
        continue
    end

    idxT = ~isnan(zT) & ~isnan(Tprof);
    idxR = ~isnan(zrho) & ~isnan(rho) & ~isnan(S) & ~isnan(Tlab);

    if sum(idxT) < 2 || sum(idxR) == 0
        cores(k).temperature_rho = T_rho;
        cores(k).brine_volume = vb;
        cores(k).gas_volume_lab = vg_pr_out;
        cores(k).gas_volume_connected = vg_connected;
        cores(k).gas_volume_disconnected = vg_disconnected;
        cores(k).density_insitu_connected = rho_si_connected;
        cores(k).density_insitu_disconnected = rho_si_disconnected;
        cores(k).temperature_bulk = NaN;
        cores(k).density_insitu_connected_bulk = NaN;
        cores(k).density_insitu_disconnected_bulk = NaN;
        continue
    end

    zTok = zT(idxT);
    Tok = min(-0.1,Tprof(idxT));

    [zTuniq,~,ic] = unique(zTok);
    Tuniq = accumarray(ic,Tok,[],@mean);

    if numel(zTuniq) < 2
        cores(k).temperature_rho = T_rho;
        cores(k).brine_volume = vb;
        cores(k).gas_volume_lab = vg_pr_out;
        cores(k).gas_volume_connected = vg_connected;
        cores(k).gas_volume_disconnected = vg_disconnected;
        cores(k).density_insitu_connected = rho_si_connected;
        cores(k).density_insitu_disconnected = rho_si_disconnected;
        cores(k).temperature_bulk = NaN;
        cores(k).density_insitu_connected_bulk = NaN;
        cores(k).density_insitu_disconnected_bulk = NaN;
        continue
    end

    T_rho(idxR) = interp1(zTuniq,Tuniq,zrho(idxR),'linear','extrap');

    TlabR = Tlab(idxR);
    Ti = T_rho(idxR);
    Si = S(idxR);
    rho_lab = rho(idxR);

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

    A = ...
        (rhoi_i ./ rhoi_pr) .* ...
        (F3_pr ./ F3_i) .* ...
        (F1_pr ./ F1_i);

    idx0 = Si == 0;
    A(idx0) = 1;

    vg_conn = max(0, ...
        1 - (1 - vg_pr) .* A);

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
    vg_pr_out(idxR) = vg_pr;
    vg_connected(idxR) = vg_conn;
    vg_disconnected(idxR) = vg_disc;
    rho_si_connected(idxR) = rho_conn;
    rho_si_disconnected(idxR) = rho_disc;

    cores(k).temperature_rho = T_rho;
    cores(k).brine_volume = vb;
    cores(k).gas_volume_lab = vg_pr_out;
    cores(k).gas_volume_connected = vg_connected;
    cores(k).gas_volume_disconnected = vg_disconnected;
    cores(k).density_insitu_connected = rho_si_connected;
    cores(k).density_insitu_disconnected = rho_si_disconnected;
    cores(k).temperature_bulk = mean(T_rho,'omitnan');
    cores(k).density_insitu_connected_bulk = mean(rho_si_connected,'omitnan');
    cores(k).density_insitu_disconnected_bulk = mean(rho_si_disconnected,'omitnan');
end

n = numel(cores);

rho_sud = nan(n,1);
rho_si_connected_sud = nan(n,1);
rho_si_disconnected_sud = nan(n,1);
hi_sud = nan(n,1);
hs_sud = nan(n,1);
Ti_sud = nan(n,1);
t_sud = NaT(n,1);
lat_sud = nan(n,1);
lon_sud = nan(n,1);
source_sud = strings(n,1);
dataset_sud = strings(n,1);
file_sud = strings(n,1);

for k = 1:n
    file_sud(k) = string(cores(k).file);
    source_sud(k) = string(cores(k).source);
    dataset_sud(k) = "SUDARCO";

    rho_sud(k) = mean(cores(k).density,'omitnan');
    rho_si_connected_sud(k) = mean(cores(k).density_insitu_connected,'omitnan');
    rho_si_disconnected_sud(k) = mean(cores(k).density_insitu_disconnected,'omitnan');

    hi_sud(k) = cores(k).ice_thickness / 100;
    hs_sud(k) = cores(k).snow_thickness / 100;
    Ti_sud(k) = mean(cores(k).temperature_rho,'omitnan');

    t_sud(k) = cores(k).date;
    lat_sud(k) = cores(k).lat;
    lon_sud(k) = cores(k).lon;
end

ice_age = repmat("SYI",n,1);

idxMOSAiC = file_sud == "MOSAiC_ice_station_14.06.2020.xlsx";
ice_age(idxMOSAiC) = "FYI";
dataset_sud(idxMOSAiC) = "MOSAiC FYI";

ice_age(file_sud == "AO2023-I-IceStation_2_day1_with_density_final.xlsx") = "FYI";
ice_age(file_sud == "AO2023-I-IceStation_3_day1_with_density_final.xlsx") = "MYI";
ice_age(file_sud == "AO2023-I-IceStation_4_day1_with_density_final.xlsx") = "MYI";
ice_age(file_sud == "AO2023-I-IceStation_5_day1_with_density_final.xlsx") = "FYI";

ice_age(file_sud == "Polhavet2022_Ice_Station_1_NP_ver110323_with_density_stratigraphy_final.xlsx") = "FYI";
ice_age(file_sud == "Polhavet2022_Ice_Station_6_Amundsen_Basin_ver150223_with_density.xlsx") = "SYI";
ice_age(file_sud == "Polhavet2022_Ice_Station_8_Nansen_Basin_ver030223_with_density_final.xlsx") = "SYI";
ice_age(file_sud == "Polhavet2022_Ice_Station_10_ver291123with_density_final.xlsx") = "FYI";

SUDARCO_core_summary = table( ...
    file_sud, ...
    dataset_sud, ...
    t_sud, ...
    ice_age, ...
    rho_sud, ...
    rho_si_connected_sud, ...
    rho_si_disconnected_sud, ...
    hi_sud, ...
    hs_sud, ...
    Ti_sud, ...
    lat_sud, ...
    lon_sud, ...
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

SUDARCO_core_summary = sortrows(SUDARCO_core_summary,{'dataset','date'});

fprintf('Imported SUDARCO dataset.\n')
fprintf('Processed %d ice cores.\n', ...
    height(SUDARCO_core_summary))
fprintf('Output saved to:\n%s\n', ...
    outputFile)

if ismember('snow_depth_m',SUDARCO_core_summary.Properties.VariableNames)
    SUDARCO_core_summary.Properties.VariableNames{'snow_depth_m'} = ...
        'snow_thickness_m';
end

if ismember('file',SUDARCO_core_summary.Properties.VariableNames) && ...
        ~ismember('core',SUDARCO_core_summary.Properties.VariableNames)
    SUDARCO_core_summary.core = SUDARCO_core_summary.file;
    SUDARCO_core_summary = movevars(SUDARCO_core_summary,'core','After','dataset');
end

Summary_SUDARCO = SUDARCO_core_summary;

save(outputFile,'Summary_SUDARCO')

%% Helper functions

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