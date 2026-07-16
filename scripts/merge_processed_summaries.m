function merge_processed_summaries(repoRoot)
%
% MERGE_PROCESSED_SUMMARIES
%
% Merges all processed campaign summary tables into the common final
% matchup table used for SnowModel-LG extraction and density analysis.
%
% The script standardizes variable names, harmonizes table schemas, applies
% selected snow-thickness corrections, and initializes placeholder columns
% for ERA5 and MERRA2 SnowModel-LG matchups.
%
% Inputs:
%   data/processed/Summary_*.mat
%
% Outputs:
%   data/final/snow_model_matchups.mat
%

close all

processedDir = fullfile(repoRoot, 'data', 'processed');
finalDir = fullfile(repoRoot, 'data', 'final');

if ~exist(finalDir, 'dir')
    mkdir(finalDir)
end

dataFolder = processedDir;

% Processed campaign summaries to merge.
summaryFiles = { ...
    fullfile(dataFolder,'Summary_MOSAIC_FYI.mat')
    fullfile(dataFolder,'Summary_MOSAIC_SYI.mat')
    fullfile(dataFolder,'Summary_MOSAiC_Leg5.mat')
    fullfile(dataFolder,'Summary_NL.mat')
    fullfile(dataFolder,'Summary_SUDARCO.mat')
    fullfile(dataFolder,'Summary_Svalbard.mat')
    fullfile(dataFolder,'Summary_CONTRASTS.mat')
    fullfile(dataFolder,'Summary_GoNorth.mat')
};

averageGoNorth = true;

% Load all available campaign summary tables.
All = table();

for i = 1:numel(summaryFiles)

    if ~isfile(summaryFiles{i})
        warning('Missing file: %s',summaryFiles{i})
        continue
    end

    S = load(summaryFiles{i});
    fn = fieldnames(S);

    Ttmp = S.(fn{1});

    % Standardize key variable types after concatenation.

    if ismember('snow_depth_m',Ttmp.Properties.VariableNames)
        Ttmp.Properties.VariableNames{'snow_depth_m'} = ...
            'snow_thickness_m';
    end

    if ismember('file',Ttmp.Properties.VariableNames) && ...
            ~ismember('core',Ttmp.Properties.VariableNames)

        Ttmp.core = string(Ttmp.file);
    end

    % variables wanted

    wantedVars = { ...
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
        'lon_degE'};

    % skip incomplete tables

    missingVars = ...
        setdiff(wantedVars,Ttmp.Properties.VariableNames);

    if ~isempty(missingVars)

        warning('Skipping %s because missing: %s', ...
            summaryFiles{i}, ...
            strjoin(missingVars,', '))

        continue
    end

    Ttmp = Ttmp(:,wantedVars);

    All = [All; Ttmp];

end

if isempty(All)
    error('No valid summary tables loaded.')
end

% STANDARDIZE TYPES

All.dataset = string(All.dataset);
All.core = string(All.core);
All.ice_age = upper(string(All.ice_age));

if ~isdatetime(All.date)
    All.date = datetime(All.date);
end

% Optionally average GoNorth observations to one FYI and one SYI record.

if averageGoNorth

    idxGN = contains(lower(All.dataset),'gonorth');

    GN = All(idxGN,:);
    All = All(~idxGN,:);

    if ~isempty(GN)

        GN.is_old = GN.ice_age ~= "FYI";

        [G,oldFlag] = findgroups(GN.is_old);

        GNavg = table();

        nG = numel(oldFlag);

        GNavg.dataset = strings(nG,1);
        GNavg.core = strings(nG,1);
        GNavg.date = NaT(nG,1);
        GNavg.ice_age = strings(nG,1);

        GNavg.temperature_C = ...
            splitapply(@(x) mean(x,'omitnan'), ...
            GN.temperature_C,G);

        GNavg.density_lab_kgm3 = ...
            splitapply(@(x) mean(x,'omitnan'), ...
            GN.density_lab_kgm3,G);

        GNavg.density_insitu_connected_kgm3 = ...
            splitapply(@(x) mean(x,'omitnan'), ...
            GN.density_insitu_connected_kgm3,G);

        GNavg.density_insitu_disconnected_kgm3 = ...
            splitapply(@(x) mean(x,'omitnan'), ...
            GN.density_insitu_disconnected_kgm3,G);

        GNavg.ice_thickness_m = ...
            splitapply(@(x) mean(x,'omitnan'), ...
            GN.ice_thickness_m,G);

        GNavg.snow_thickness_m = ...
            splitapply(@(x) mean(x,'omitnan'), ...
            GN.snow_thickness_m,G);

        GNavg.lat_degN = ...
            splitapply(@(x) mean(x,'omitnan'), ...
            GN.lat_degN,G);

        GNavg.lon_degE = ...
            splitapply(@(x) mean(x,'omitnan'), ...
            GN.lon_degE,G);

        for k = 1:nG

            idxGroup = G == k;

            GNavg.date(k) = ...
                mean(GN.date(idxGroup),'omitnan');

            if oldFlag(k)

                GNavg.dataset(k) = "GoNorth SYI";
                GNavg.core(k) = "GoNorth SYI";
                GNavg.ice_age(k) = "SYI";

            else

                GNavg.dataset(k) = "GoNorth FYI";
                GNavg.core(k) = "GoNorth FYI";
                GNavg.ice_age(k) = "FYI";

            end
        end

        All = [All; GNavg];

    end
end

% Build final common matchup table.

D = table();

D.dataset = All.dataset;
D.core = All.core;
D.date = All.date;

D.lat = double(All.lat_degN);
D.lon = double(All.lon_degE);

D.hs = double(All.snow_thickness_m);

D.ice_age = string(All.ice_age);

D.temperature_C = double(All.temperature_C);

D.ice_thickness_m = ...
    double(All.ice_thickness_m);

D.density_lab_kgm3 = ...
    double(All.density_lab_kgm3);

D.density_insitu_connected_kgm3 = ...
    double(All.density_insitu_connected_kgm3);

D.density_insitu_disconnected_kgm3 = ...
    double(All.density_insitu_disconnected_kgm3);

% Remove rows without valid date or geographic coordinates.

validD = ...
    ~isnat(D.date) & ...
    ~isnan(D.lat) & ...
    ~isnan(D.lon);

D = D(validD,:);

% Apply dataset-specific snow-thickness corrections.

md = month(D.date);
dd = day(D.date);

% Set summer snow thickness to zero for late June through August
% observations where snow cover is assumed absent.
summer_mask = ...
    (md == 6 & dd >= 25) | ...
    md == 7 | ...
    (md == 8 & dd <= 31);

D.hs(summer_mask) = 0;

% Fill selected missing MOSAiC SYI snow-thickness values using the mean of
% nearby available MOSAiC SYI observations.
idx_fill = [31 33 34];

idx_fill = idx_fill(idx_fill <= height(D));

if ~isempty(idx_fill) && height(D) >= 30

    hs_ref = mean(D.hs(22:30),'omitnan');

    D.hs(idx_fill) = hs_ref;

end

% Initialize SnowModel-LG matchup variables. These columns are populated
% later by add_SMLG_model_snow.m after ERA5- and MERRA2-forced SnowModel-LG
% snow depths have been matched to each observation.

D.hs_smlg = nan(height(D),1);
D.t_smlg = NaT(height(D),1);

D.hs_merra = nan(height(D),1);
D.t_merra = NaT(height(D),1);

% SAVE

outputFile = fullfile(finalDir, ...
    'snow_model_matchups.mat');

save(outputFile,'D','All')

fprintf('Merged processed summaries.\n')
fprintf('Saved %d rows to:\n%s\n', ...
    height(D), outputFile)

end