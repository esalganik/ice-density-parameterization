clear; clc

%% Repository setup

% Repository root directory
repoRoot = fileparts(mfilename('fullpath'));

% Add repository scripts and external packages to MATLAB path
addpath(genpath(fullfile(repoRoot, 'scripts')))
addpath(genpath(fullfile(repoRoot, 'external')))

%% Process observational datasets

CONTRASTS_processing(repoRoot)      % CONTRASTS
GoNorth_processing(repoRoot)        % GoNorth
MOSAIC_FYI_processing(repoRoot)     % MOSAiC FYI
MOSAiC_leg5_processing(repoRoot)    % MOSAiC Leg 5
MOSAiC_SYI_processing(repoRoot)     % MOSAiC SYI
Nansen_Legacy_processing(repoRoot)  % Nansen Legacy
SUDARCO_processing(repoRoot)        % SUDARCO
Svalbard_processing(repoRoot)       % Svalbard

%% Merge processed datasets

merge_processed_summaries(repoRoot)

%% Add external snow products

finalDir = fullfile(repoRoot, 'data', 'final');

modelFile = fullfile(finalDir, ...
    'snow_model_matchups_with_models.mat');

era5File = fullfile(repoRoot, 'data', 'model', 'ERA5', ...
    'SM_snod_ERA5_01Aug1980-31Jul2021_v01.nc');

merraFiles = { ...
    fullfile(repoRoot, 'data', 'model', 'MERRA2', ...
    'SM_snod_MERRA2_ease_01Aug1980-31Jul2018.nc')
    fullfile(repoRoot, 'data', 'model', 'MERRA2', ...
    'SM_snod_MERRA2_ease_01Aug2018-31Jul2021.nc')
    fullfile(repoRoot, 'data', 'model', 'MERRA2', ...
    'SM_snod_MERRA2_ease_01Aug2021-31Jul2022.nc')
    fullfile(repoRoot, 'data', 'model', 'MERRA2', ...
    'SM_snod_MERRA2_ease_01Aug2022-31Jul2023.nc')
    };

requiredModelFiles = [{era5File}; merraFiles];

haveModelFiles = all(cellfun(@(f) exist(f,'file'), ...
    requiredModelFiles));

if haveModelFiles

    fprintf('SM-LG source files found. Regenerating snow matchups.\n')
    add_SMLG_model_snow(repoRoot)

elseif exist(modelFile,'file')

    fprintf(['SM-LG source files missing. Using existing ' ...
        'processed matchup file:\n%s\n'], modelFile)

else

    error(['SM-LG source files are missing and no processed ' ...
        'matchup file was found.\n' ...
        'Provide the NetCDF files in data/model/ or provide:\n%s'], ...
        modelFile)

end

%% Generate manuscript figures

% Reuse existing map if available to avoid requiring the large ETOPO1 topography dataset for reruns.
mapFile = fullfile(repoRoot, ...
    'figures', 'Map_density_lab.png');

if exist(mapFile,'file')
    fprintf('Using existing density map:\n%s\n', mapFile)
else
    plot_density_map(repoRoot)
end

analyze_density_model(repoRoot)     % Figure 1
analyze_snow_dependence(repoRoot)   % Figure 2