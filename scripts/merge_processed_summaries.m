function merge_processed_summaries(repoRoot)

close all

processedDir = fullfile(repoRoot, 'data', 'processed');
finalDir = fullfile(repoRoot, 'data', 'final');

if ~exist(finalDir, 'dir')
    mkdir(finalDir)
end

dataFolder = processedDir;

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

% LOAD ALL TABLES

All = table();

for i = 1:numel(summaryFiles)

    if ~isfile(summaryFiles{i})
        warning('Missing file: %s',summaryFiles{i})
        continue
    end

    S = load(summaryFiles{i});
    fn = fieldnames(S);

    Ttmp = S.(fn{1});

    % standardize names

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

% AVERAGE GONORTH

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

% BUILD FINAL TABLE

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

% REMOVE BAD GEO/DATE

validD = ...
    ~isnat(D.date) & ...
    ~isnan(D.lat) & ...
    ~isnan(D.lon);

D = D(validD,:);

% SNOW THICKNESS FIXES

md = month(D.date);
dd = day(D.date);

summer_mask = ...
    (md == 6 & dd >= 25) | ...
    md == 7 | ...
    (md == 8 & dd <= 31);

D.hs(summer_mask) = 0;

% fill missing MOSAiC SYI snow

idx_fill = [31 33 34];

idx_fill = idx_fill(idx_fill <= height(D));

if ~isempty(idx_fill) && height(D) >= 30

    hs_ref = mean(D.hs(22:30),'omitnan');

    D.hs(idx_fill) = hs_ref;

end

% placeholders

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