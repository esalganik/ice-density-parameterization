function make_mosaic_coring_density_figure(repoRoot)

close all

rawDir = fullfile(repoRoot, 'data', 'raw');
figureDir = fullfile(repoRoot, 'figures');

dataFile = fullfile(rawDir, 'MOSAiC_coring_freeboards.xlsx');

rho_w = 1024;

rhoFilterMin = 650;
rhoFilterMax = 960;

rhoPlotMin = 840;
rhoPlotMax = 960;

sheets = {'FYI','SYI'};
iceTypes = {'FYI','SYI'};

results = table();

for s = 1:numel(sheets)

    C = readcell(dataFile, 'Sheet', sheets{s});

    nCols = size(C,2);

    thickCols = 1:2:nCols;
    draftCols = 2:2:nCols;

    nPairs = min(numel(thickCols), numel(draftCols));

    for k = 1:nPairs

        cH = thickCols(k);
        cD = draftCols(k);

        coreDate = C{2,cH};
        snowH = C{5,cH};

        if isempty(coreDate) || isempty(snowH)
            continue
        end

        coreDate = toDatetime(coreDate);
        snowH = toDouble(snowH);

        iceH = cellfun(@toDouble, C(8:end,cH));
        iceD = cellfun(@toDouble, C(8:end,cD));

        snowH = snowH ./ 100;
        iceH  = iceH ./ 100;
        iceD  = iceD ./ 100;

        good0 = ~isnan(iceH) & ~isnan(iceD) & iceH > 0;

        tmp = table();

        tmp.CoreDate = repmat(coreDate, sum(good0), 1);
        tmp.IceType  = repmat(string(iceTypes{s}), sum(good0), 1);

        tmp.SnowThickness = repmat(snowH, sum(good0), 1);
        tmp.IceThickness  = iceH(good0);
        tmp.IceDraft      = iceD(good0);
        tmp.Freeboard     = iceH(good0) - iceD(good0);

        results = [results; tmp];

    end
end

yearVec = year(results.CoreDate);
aug1 = datetime(yearVec, 8, 1);
tA1 = days(results.CoreDate - aug1);

results.SnowDensity = 0.35 .* tA1 + 239.78;

results.IceDensityRaw = ...
    rho_w .* (results.IceDraft ./ results.IceThickness) ...
    - results.SnowDensity .* ...
    (results.SnowThickness ./ results.IceThickness);

results.IceDensity = results.IceDensityRaw;

bad = results.IceDensity < rhoFilterMin | ...
      results.IceDensity > rhoFilterMax;

results.IceDensity(bad) = NaN;

good = ~isnan(results.IceDensity);

isFYI = results.IceType == "FYI";
isSYI = results.IceType == "SYI";

isBeforeJuly = month(results.CoreDate) < 7;
isJulyOnward = month(results.CoreDate) >= 7;

rhoAll = results.IceDensity(good);
rhoFYI = results.IceDensity(good & isFYI);
rhoSYI = results.IceDensity(good & isSYI);

meanAll = mean(rhoAll, 'omitnan');
meanFYI = mean(rhoFYI, 'omitnan');
meanSYI = mean(rhoSYI, 'omitnan');

meanFYIBeforeJuly = mean(results.IceDensity(good & isFYI & isBeforeJuly), 'omitnan');
meanFYIJulyOnward = mean(results.IceDensity(good & isFYI & isJulyOnward), 'omitnan');

meanSYIBeforeJuly = mean(results.IceDensity(good & isSYI & isBeforeJuly), 'omitnan');
meanSYIJulyOnward = mean(results.IceDensity(good & isSYI & isJulyOnward), 'omitnan');

nFYIBeforeJuly = sum(good & isFYI & isBeforeJuly);
nFYIJulyOnward = sum(good & isFYI & isJulyOnward);

nSYIBeforeJuly = sum(good & isSYI & isBeforeJuly);
nSYIJulyOnward = sum(good & isSYI & isJulyOnward);

fig = figure;
set(fig, 'Units', 'inches', 'Position', [0.5 5 10 4])

tl = tiledlayout(fig, 1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax1 = nexttile(tl, 1);
hold(ax1, 'on')

scatter(ax1, results.CoreDate(good & isFYI), ...
        results.IceDensity(good & isFYI), ...
        30, ...
        'filled', ...
        'MarkerFaceColor', [0 0.4470 0.7410], ...
        'MarkerFaceAlpha', 0.25, ...
        'MarkerEdgeAlpha', 0.25)

scatter(ax1, results.CoreDate(good & isSYI), ...
        results.IceDensity(good & isSYI), ...
        30, ...
        'filled', ...
        'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
        'MarkerFaceAlpha', 0.25, ...
        'MarkerEdgeAlpha', 0.25)

dateFYI = dateshift(results.CoreDate(good & isFYI), 'start', 'day');
dateSYI = dateshift(results.CoreDate(good & isSYI), 'start', 'day');

[GFYI, dFYI] = findgroups(dateFYI);
[GSYI, dSYI] = findgroups(dateSYI);

meanDateFYI = splitapply(@mean, results.IceDensity(good & isFYI), GFYI);
meanDateSYI = splitapply(@mean, results.IceDensity(good & isSYI), GSYI);

plot(ax1, dFYI, meanDateFYI, ...
    'o', ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', [0 0.4470 0.7410], ...
    'LineWidth', 1.8)

plot(ax1, dSYI, meanDateSYI, ...
    'o', ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', [0.8500 0.3250 0.0980], ...
    'LineWidth', 1.8)

xline(ax1, datetime(2020,7,1), ':', ...
    'Color', [0.4 0.4 0.4], ...
    'LineWidth', 1.2)

xPre1 = datetime(2019,10,1);
xPre2 = datetime(2020,7,1);
xJul1 = datetime(2020,7,1);
xJul2 = datetime(2020,8,15);

plot(ax1, [xPre1 xPre2], [meanFYIBeforeJuly meanFYIBeforeJuly], ...
    '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.8)

plot(ax1, [xJul1 xJul2], [meanFYIJulyOnward meanFYIJulyOnward], ...
    '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.8)

plot(ax1, [xPre1 xPre2], [meanSYIBeforeJuly meanSYIBeforeJuly], ...
    '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.8)

plot(ax1, [xJul1 xJul2], [meanSYIJulyOnward meanSYIJulyOnward], ...
    '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.8)

txtPreJuly = sprintf([ ...
    'pre-July mean:\n' ...
    'FYI %.0f, N=%d\n' ...
    'SYI %.0f, N=%d'], ...
    meanFYIBeforeJuly, nFYIBeforeJuly, ...
    meanSYIBeforeJuly, nSYIBeforeJuly);

txtJuly = sprintf([ ...
    'July mean:\n' ...
    'FYI %.0f, N=%d\n' ...
    'SYI %.0f, N=%d'], ...
    meanFYIJulyOnward, nFYIJulyOnward, ...
    meanSYIJulyOnward, nSYIJulyOnward);

text(ax1, datetime(2020,2,1), 858, txtPreJuly, ...
    'FontSize', 8, ...
    'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', ...
    'Margin', 2)

text(ax1, datetime(2020,4,20), 858, txtJuly, ...
    'FontSize', 8, ...
    'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', ...
    'Margin', 2)

ylim(ax1, [rhoPlotMin rhoPlotMax])
xlim(ax1, [datetime(2019,10,1) datetime(2020,8,15)])

ylabel(ax1, 'Effective sea-ice density (kg m^{-3})')
xlabel(ax1, 'Date')
title(ax1, 'Hydrostatically derived effective sea-ice density', ...
    'FontWeight', 'normal')

xtickformat(ax1, 'MMM')
box(ax1, 'on')

lgd1 = legend(ax1,...
 {'FYI cores',...
  'SYI cores',...
  'FYI event mean',...
  'SYI event mean',...
  'Seasonal means'},...
 'Location','southwest');

lgd1.ItemTokenSize = [14 10];

ax2 = nexttile(tl, 2);
hold(ax2, 'on')

binWidth = 2 * iqr(rhoAll) / nthroot(numel(rhoAll), 3);
binWidth = round(binWidth / 10) * 10;

if binWidth <= 0
    binWidth = 10;
end

edges = rhoPlotMin:binWidth:rhoPlotMax;

histogram(ax2, rhoFYI, edges, ...
    'FaceColor', [0 0.4470 0.7410], ...
    'FaceAlpha', 0.35)

histogram(ax2, rhoSYI, edges, ...
    'FaceColor', [0.8500 0.3250 0.0980], ...
    'FaceAlpha', 0.35)

[xiFYI, fFYI] = getKDE(rhoFYI);
[xiSYI, fSYI] = getKDE(rhoSYI);

plot(ax2, xiFYI, fFYI * binWidth * numel(rhoFYI), ...
    'Color', [0 0.4470 0.7410], ...
    'LineWidth', 2)

plot(ax2, xiSYI, fSYI * binWidth * numel(rhoSYI), ...
    'Color', [0.8500 0.3250 0.0980], ...
    'LineWidth', 2)

xline(ax2, meanFYI, '--', ...
    'Color', [0 0.4470 0.7410], ...
    'LineWidth', 2)

xline(ax2, meanSYI, '--', ...
    'Color', [0.8500 0.3250 0.0980], ...
    'LineWidth', 2)

xlabel(ax2, 'Effective sea-ice density (kg m^{-3})')
ylabel(ax2, 'Counts')

title(ax2, 'Distribution of effective sea-ice density', ...
    'FontWeight', 'normal')

xlim(ax2, [rhoPlotMin rhoPlotMax])

yl = ylim(ax2);
ylim(ax2, [0 ceil(yl(2)/10)*10 + 10])

legend(ax2, ...
    {'FYI', ...
     'SYI', ...
     'FYI KDE', ...
     'SYI KDE', ...
     'FYI mean', ...
     'SYI mean'}, ...
    'Location', 'northwest', ...
    'Box', 'off')

box(ax2, 'on')

fsz = 11;

text(ax1, -0.13, 1.035, '(a)', ...
    'Units', 'normalized', ...
    'FontSize', fsz-1)

text(ax2, -0.13, 1.035, '(b)', ...
    'Units', 'normalized', ...
    'FontSize', fsz-1)

set(findall(fig, 'Type', 'axes'), 'Toolbar', [])

outputFigure = fullfile(figureDir, 'FigS1.png');
exportgraphics(fig, outputFigure, 'Resolution', 300)

close(fig)

fprintf('Generated MOSAiC coring density figure.\n')
fprintf('Upper density filter: %.0f kg m-3\n', rhoFilterMax)
fprintf('Saved figure to:\n%s\n', outputFigure)

end

function x = toDouble(v)

if isnumeric(v)
    x = double(v);

elseif ischar(v) || isstring(v)
    x = str2double(strrep(string(v), ',', '.'));

else
    x = NaN;
end

end

function d = toDatetime(v)

if isdatetime(v)

    d = v;

elseif isnumeric(v)

    d = datetime(v, 'ConvertFrom', 'excel');

elseif ischar(v) || isstring(v)

    d = datetime(v);

else

    d = NaT;

end

end

function [xi, f] = getKDE(x)

x = x(~isnan(x));

if numel(x) < 2
    xi = NaN;
    f = NaN;
    return
end

[f, xi] = ksdensity(x);

end