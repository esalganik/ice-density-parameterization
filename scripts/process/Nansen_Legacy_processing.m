function Nansen_Legacy_processing(repoRoot)

close all

rawDir = fullfile(repoRoot, 'data', 'raw');
processedDir = fullfile(repoRoot, 'data', 'processed');

if ~exist(processedDir, 'dir')
    mkdir(processedDir)
end

inputFile = fullfile(rawDir, ...
    'NL_density.xlsx');

outputFile = fullfile(processedDir, ...
    'Summary_NL.mat');

T1 = readtable(inputFile,'Sheet','density');
T2 = readtable(inputFile,'Sheet','temperature');

n = table2array(T1(:,1));
rho_all = table2array(T1(:,19));
t_all = table2array(T1(:,2));
zzrho_all = table2array(T1(:,9:10));
T_lab_all = table2array(T1(:,11));
Srho_all = table2array(T1(:,12));
hrho_all = table2array(T1(:,8));
drho_all = table2array(T1(:,8)) - table2array(T1(:,6));
vg_orig_all = table2array(T1(:,14));
lon_all = table2array(T1(:,4));
lat_all = table2array(T1(:,5));
sn_all = table2array(T1(:,7));

n2 = table2array(T2(:,1));
hT_all = table2array(T2(:,8));
dT_all = table2array(T2(:,8)) - table2array(T2(:,6));
zT_all = table2array(T2(:,9));
T_all = table2array(T2(:,11));

for i = 1:max(n)
    rho{i} = rho_all(n == i);
    zzrho{i} = zzrho_all(n == i,:);
    Srho{i} = Srho_all(n == i);
    T_lab{i} = T_lab_all(n == i)*0 - 20;

    t0{i} = t_all(n == i);
    t(i) = t0{i}(1);

    drho0{i} = drho_all(n == i);
    drho(i) = drho0{i}(1);

    vg_orig{i} = vg_orig_all(n == i);

    hrho0{i} = hrho_all(n == i);
    hrho(i) = hrho0{i}(1);

    snrho0{i} = sn_all(n == i);
    snrho(i) = snrho0{i}(1);

    lon0{i} = lon_all(n == i);
    lat0{i} = lat_all(n == i);
    lon(i) = lon0{i}(1);
    lat(i) = lat0{i}(1);

    T{i} = T_all(n2 == i);
    zT{i} = zT_all(n2 == i);
    dT{i} = dT_all(n2 == i);

    hT0{i} = hT_all(n2 == i);
    hT(i) = hT0{i}(1);
end

for i = 24:27
    t0{i} = t0{i} - 365;
    t(i) = t(i) - 365;
end

for i = 1:max(n)

    zrho{i} = mean(zzrho{i},2) * hrho(i) / zzrho{i}(end,2);

    T_profile = min(-0.1,T{i});
    zT_scaled = zT{i} * hT(i) / zT{i}(end);

    okT = ~isnan(zT_scaled) & ~isnan(T_profile);

    [zT_unique,~,ic] = unique(zT_scaled(okT));
    T_unique = accumarray(ic,T_profile(okT),[],@mean);

    T_rho{i} = interp1(zT_unique,T_unique,zrho{i},'linear','extrap');

    Tlab = T_lab{i};
    Ti = T_rho{i};
    Si = Srho{i};
    rho_lab = rho{i};

    [F1_pr,F2_pr] = F1F2_seaice(Tlab);
    [F1_i,F2_i] = F1F2_seaice(Ti);

    rhoi_pr = 917 - 0.1403*Tlab;
    rhoi_i = 917 - 0.1403*Ti;

    vb_pr = rho_lab .* Si ./ F1_pr;

    vg_pr{i} = max(0, ...
        1 - rho_lab .* ...
        (F1_pr - rhoi_pr .* Si/1000 .* F2_pr) ./ ...
        (rhoi_pr .* F1_pr));

    F3_pr = rhoi_pr .* Si/1000 ./ ...
        (F1_pr - rhoi_pr .* Si/1000 .* F2_pr);

    F3_i = rhoi_i .* Si/1000 ./ ...
        (F1_i - rhoi_i .* Si/1000 .* F2_i);

    vb_raw = vb_pr .* F1_pr ./ F1_i / 1000;

    vb_rho{i} = vb_raw;
    vb_rho{i}(vb_rho{i} > 0.6 | vb_rho{i} < 0) = NaN;

    vb_calc = vb_raw;
    vb_calc(vb_calc > 0.6 | vb_calc < 0) = 0.6;

    A = ...
        (rhoi_i ./ rhoi_pr) .* ...
        (F3_pr ./ F3_i) .* ...
        (F1_pr ./ F1_i);

    idx0 = Si == 0;
    A(idx0) = 1;

    vg_connected{i} = max(0, ...
        1 - (1 - vg_pr{i}) .* A);

    vg_disconnected{i} = vg_pr{i};

    rho_si_connected{i} = ...
        (1 - vg_connected{i}) .* rhoi_i .* F1_i ./ ...
        (F1_i - rhoi_i .* Si/1000 .* F2_i);

    rho_si_disconnected{i} = ...
        (1 - vg_disconnected{i}) .* rhoi_i .* F1_i ./ ...
        (F1_i - rhoi_i .* Si/1000 .* F2_i);

    rho_si_connected{i}(isnan(vb_rho{i})) = NaN;
    rho_si_disconnected{i}(isnan(vb_rho{i})) = NaN;

    rho_si_bulk_connected(i) = mean(rho_si_connected{i},'omitnan');
    rho_si_bulk_disconnected(i) = mean(rho_si_disconnected{i},'omitnan');
    T_bulk(i) = mean(T_rho{i},'omitnan');
end

rho_bulk = nan(max(n),1);
rho_si_connected_bulk = nan(max(n),1);
rho_si_disconnected_bulk = nan(max(n),1);
T_bulk = nan(max(n),1);
hi = nan(max(n),1);
hs = nan(max(n),1);
date_core = NaT(max(n),1);
lat_core = nan(max(n),1);
lon_core = nan(max(n),1);

for i = 1:max(n)
    rho_bulk(i) = mean(rho{i},'omitnan');
    rho_si_connected_bulk(i) = mean(rho_si_connected{i},'omitnan');
    rho_si_disconnected_bulk(i) = mean(rho_si_disconnected{i},'omitnan');
    T_bulk(i) = mean(T_rho{i},'omitnan');

    hi(i) = hrho(i);
    hs(i) = snrho(i);
    date_core(i) = t(i);
    lat_core(i) = lat(i);
    lon_core(i) = lon(i);
end

ice_age = repmat("FYI",max(n),1);
dataset = repmat("Nansen Legacy FYI",max(n),1);

NL_core_summary = table( ...
    dataset, ...
    string((1:max(n))'), ...
    date_core, ...
    ice_age, ...
    T_bulk, ...
    rho_bulk, ...
    rho_si_connected_bulk, ...
    rho_si_disconnected_bulk, ...
    hi, ...
    hs, ...
    lat_core, ...
    lon_core, ...
    'VariableNames', { ...
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

fprintf('Imported Nansen Legacy dataset.\n')
fprintf('Processed %d ice cores.\n', ...
    height(NL_core_summary))
fprintf('Output saved to:\n%s\n', ...
    outputFile)

save(outputFile,'NL_core_summary')

%% Helpers

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