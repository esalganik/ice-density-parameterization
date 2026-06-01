clear; clc

%% Workflow switches

doProcessing = false;     % Reprocess raw observational datasets
doMerge = false;          % Merge processed datasets
doSnowMatchups = false;   % Add ERA5/MERRA2 SnowModel-LG snow products
doStatistics = true;     % Generate manuscript statistics
doFigures = false;        % Reproduce manuscript figures
doSupportingInfo = false; % Generate Supporting Information figures and tables

%% Repository setup

repoRoot = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(repoRoot, 'scripts')))
addpath(genpath(fullfile(repoRoot, 'external')))

%% Process observational datasets

if doProcessing

    CONTRASTS_processing(repoRoot)
    GoNorth_processing(repoRoot)
    MOSAIC_FYI_processing(repoRoot)
    MOSAiC_leg5_processing(repoRoot)
    MOSAiC_SYI_processing(repoRoot)
    Nansen_Legacy_processing(repoRoot)
    SUDARCO_processing(repoRoot)
    Svalbard_processing(repoRoot)

end

%% Process Supporting Information datasets

if doSupportingInfo

    Ridges_processing(repoRoot)

end

%% Merge processed datasets

if doMerge

    merge_processed_summaries(repoRoot)

end

%% Add external snow products

if doSnowMatchups

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

end

%% Generate manuscript statistics

if doStatistics

    resultsDir = fullfile(repoRoot, 'results');

    if ~exist(resultsDir, 'dir')
        mkdir(resultsDir)
    end

    statisticsFile = fullfile(resultsDir, 'manuscript_statistics.txt');

    if exist(statisticsFile,'file')
        delete(statisticsFile)
    end

    diary(statisticsFile)

    try
        manuscript_statistics(repoRoot)
    catch ME
        diary off
        rethrow(ME)
    end

    diary off

end

%% Generate manuscript figures

if doFigures

    mapFile = fullfile(repoRoot, ...
        'figures', 'Map_density_lab.png');

    if exist(mapFile,'file')
        fprintf('Using existing density map:\n%s\n', mapFile)
    else
        plot_density_map(repoRoot)
    end

    analyze_density_model(repoRoot)
    analyze_snow_dependence(repoRoot)

end

%% Generate Supporting Information figures and tables

if doSupportingInfo

    make_cross_validation_table(repoRoot)
    make_mosaic_coring_density_figure(repoRoot)
    make_si_ridge_density_figure(repoRoot)
    make_figureS3_snow_dependence_merra(repoRoot)

end