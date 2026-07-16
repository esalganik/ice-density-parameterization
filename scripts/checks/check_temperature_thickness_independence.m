%% check_temperature_thickness_independence.m
%
% Standalone diagnostic for the independent contribution of ice thickness
% to the temperature-thickness sea-ice density parameterization.
%
% The script:
%   1. Loads the same merged dataset used by analyze_density_model.m.
%   2. Quantifies the observed temperature-thickness association.
%   3. Compares a temperature-only model with the full temperature-thickness
%      model using leave-one-dataset-out (LODO) cross-validation.
%   4. Tests whether residual density variability from the temperature-only
%      model still depends on thickness, separately in warm and cold regimes.
%   5. Produces a four-panel diagnostic figure and prints manuscript-ready
%      summary statistics.
%
% Important interpretation:
% A global correlation between temperature and thickness mixes seasonal
% evolution with differences among ice types and campaigns. The main test is
% therefore the LODO comparison: does thickness improve prediction after the
% temperature transition is already represented?

clear
close all
clc

%% User settings
repoRoot = 'C:\Users\evsalg001\Documents\MATLAB\Density parametrization\script';

densityMode = 1;
% 1 = laboratory measurements
% 2 = in-situ connected porosity
% 3 = in-situ disconnected porosity

excludeSvalbard = false;   % Keep consistent with the main temperature model
saveFigure = false;

%% Paths and input
finalDir = fullfile(repoRoot,'data','final');
figureDir = fullfile(repoRoot,'figures');
inputFile = fullfile(finalDir,'snow_model_matchups_with_models.mat');

if ~exist(inputFile,'file')
    error('Input file not found:\n%s',inputFile)
end

if saveFigure && ~exist(figureDir,'dir')
    mkdir(figureDir)
end

load(inputFile,'D')

rho_cols = {"density_lab_kgm3", ...
            "density_insitu_connected_kgm3", ...
            "density_insitu_disconnected_kgm3"};
rho_names = {"laboratory","connected","disconnected"};
rho_col = rho_cols{densityMode};

requiredVars = ["date","ice_age","dataset","temperature_C", ...
                "ice_thickness_m",rho_col];
missingVars = requiredVars(~ismember(requiredVars,string(D.Properties.VariableNames)));
if ~isempty(missingVars)
    error('D is missing required variables: %s', ...
        strjoin(cellstr(missingVars),', '))
end

%% Use the same observations as the main temperature model
if ~isdatetime(D.date)
    D.date = datetime(D.date);
end

date = D.date;
T = D.temperature_C;
h = D.ice_thickness_m;
rho = D.(rho_col);
dataset_id = string(D.dataset);
ice_type = upper(string(D.ice_age));
is_svalbard = contains(lower(dataset_id),'svalbard');

valid = ~isnat(date) & isfinite(T) & isfinite(h) & isfinite(rho) & ...
        T < 0 & h > 0;
if excludeSvalbard
    valid = valid & ~is_svalbard;
end

date = date(valid);
T = T(valid);
h = h(valid);
rho = rho(valid);
dataset_id = dataset_id(valid);
ice_type = ice_type(valid);
is_svalbard = is_svalbard(valid);
N = numel(rho);

%% Model settings copied from analyze_density_model.m
h_cold_peak = 1.10;
rho_cold_peak = 911;
h_cold_high = 3.50;
rho_cold_high = 908.5;
rho_cold_0_grid = 905:1:912;

Tcrit_grid = -3.0:0.05:-2.0;
hcrit_warm_grid = 1.8:0.1:2.3;
rho_warm_plateau_grid = 895:1:905;

min_n_cold = 8;
min_n_warm = 8;

%% 1. Descriptive temperature-thickness association
[rPearson,pPearson] = corr(T,h,'Type','Pearson','Rows','complete');
[rSpearman,pSpearman] = corr(T,h,'Type','Spearman','Rows','complete');

% Remove campaign means before calculating a within-dataset association.
% This reduces differences among campaigns, although it does not fully remove
% the seasonal dependence.
T_within = nan(size(T));
h_within = nan(size(h));
groups = unique(dataset_id);
for i = 1:numel(groups)
    idx = dataset_id == groups(i);
    T_within(idx) = T(idx) - mean(T(idx),'omitnan');
    h_within(idx) = h(idx) - mean(h(idx),'omitnan');
end
[rWithin,pWithin] = corr(T_within,h_within, ...
    'Type','Pearson','Rows','complete');

%% 2. Fit full and temperature-only models to all observations
fitFull = fit_full_model(rho,T,h, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high, ...
    rho_cold_0_grid,Tcrit_grid,rho_warm_plateau_grid, ...
    hcrit_warm_grid,min_n_cold,min_n_warm);

pFull = [fitFull.Tcrit fitFull.rho_warm_plateau ...
         fitFull.warm_slope fitFull.hcrit_warm fitFull.rho_cold_0];

fitTemp = fit_temperature_only_model(rho,T,Tcrit_grid, ...
    min_n_cold,min_n_warm);
pTemp = [fitTemp.Tcrit fitTemp.rhoCold fitTemp.rhoWarm];

rhohatFull = predict_full_model(T,h,pFull, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
rhohatTemp = predict_temperature_only_model(T,pTemp);

metFull = regression_metrics(rho,rhohatFull,5);
metTemp = regression_metrics(rho,rhohatTemp,3);

%% 3. Leave-one-dataset-out cross-validation
cvFull = lodo_full_model(rho,T,h,dataset_id,pFull, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high, ...
    rho_cold_0_grid,Tcrit_grid,rho_warm_plateau_grid, ...
    hcrit_warm_grid,min_n_cold,min_n_warm);

cvTemp = lodo_temperature_only_model(rho,T,dataset_id,pTemp, ...
    Tcrit_grid,min_n_cold,min_n_warm);

common = isfinite(cvFull.yhat) & isfinite(cvTemp.yhat);
cvFullCommon = regression_metrics(rho(common),cvFull.yhat(common),5);
cvTempCommon = regression_metrics(rho(common),cvTemp.yhat(common),3);

deltaRMSE = cvTempCommon.rmse - cvFullCommon.rmse;
deltaR2 = cvFullCommon.r2 - cvTempCommon.r2;

%% 4. Does thickness explain temperature-only residuals?
% Use cross-validated residuals to avoid evaluating the temperature-only
% model on the same observations used for fitting.
residTempCV = rho - cvTemp.yhat;
TcritReference = fitFull.Tcrit;
isWarm = T > TcritReference & isfinite(residTempCV);
isCold = T <= TcritReference & isfinite(residTempCV);

warmStats = residual_thickness_test(h(isWarm),residTempCV(isWarm));
coldStats = residual_thickness_test(h(isCold),residTempCV(isCold));
allStats = residual_thickness_test(h(common),residTempCV(common));

%% 5. Fold-level model comparison
foldTable = build_fold_table(rho,dataset_id,cvTemp.yhat,cvFull.yhat);
validFolds = isfinite(foldTable.RMSE_temperature_only) & ...
             isfinite(foldTable.RMSE_full);

if sum(validFolds) >= 2
    [pFold,~,foldTestStats] = signrank( ...
        foldTable.RMSE_temperature_only(validFolds), ...
        foldTable.RMSE_full(validFolds));
else
    pFold = NaN;
    foldTestStats = struct('signedrank',NaN);
end

%% Printed results
fprintf('\n============================================================\n')
fprintf('TEMPERATURE-THICKNESS INDEPENDENCE DIAGNOSTIC\n')
fprintf('Density definition: %s\n',rho_names{densityMode})
fprintf('N = %d; datasets = %d; Svalbard excluded = %d\n', ...
    N,numel(groups),excludeSvalbard)
fprintf('============================================================\n')

fprintf('\nObserved association between temperature and thickness:\n')
fprintf('  Global Pearson r    = %+.3f, p = %.3g\n',rPearson,pPearson)
fprintf('  Global Spearman rho = %+.3f, p = %.3g\n',rSpearman,pSpearman)
fprintf('  Within-dataset r    = %+.3f, p = %.3g\n',rWithin,pWithin)

fprintf('\nIn-sample model comparison:\n')
fprintf('  Temperature only:       R2 = %.3f, RMSE = %.2f kg m-3\n', ...
    metTemp.r2,metTemp.rmse)
fprintf('  Temperature + thickness: R2 = %.3f, RMSE = %.2f kg m-3\n', ...
    metFull.r2,metFull.rmse)

fprintf('\nLODO cross-validation on common predictions (N = %d):\n',sum(common))
fprintf('  Temperature only:       R2 = %.3f, RMSE = %.2f kg m-3\n', ...
    cvTempCommon.r2,cvTempCommon.rmse)
fprintf('  Temperature + thickness: R2 = %.3f, RMSE = %.2f kg m-3\n', ...
    cvFullCommon.r2,cvFullCommon.rmse)
fprintf('  Improvement from thickness: Delta RMSE = %.2f kg m-3; Delta R2 = %.3f\n', ...
    deltaRMSE,deltaR2)

fprintf('\nThickness dependence of cross-validated temperature-only residuals:\n')
fprintf('  All data: slope = %+.2f kg m-3 m-1, r = %+.3f, p = %.3g, N = %d\n', ...
    allStats.slope,allStats.r,allStats.p,allStats.N)
fprintf('  Warm regime (T > %.2f C): slope = %+.2f kg m-3 m-1, r = %+.3f, p = %.3g, N = %d\n', ...
    TcritReference,warmStats.slope,warmStats.r,warmStats.p,warmStats.N)
fprintf('  Cold regime (T <= %.2f C): slope = %+.2f kg m-3 m-1, r = %+.3f, p = %.3g, N = %d\n', ...
    TcritReference,coldStats.slope,coldStats.r,coldStats.p,coldStats.N)

fprintf('\nFold-level paired comparison:\n')
fprintf('  Median RMSE improvement = %.2f kg m-3\n', ...
    median(foldTable.RMSE_temperature_only(validFolds) - ...
           foldTable.RMSE_full(validFolds),'omitnan'))
fprintf('  Wilcoxon signed-rank p = %.3g\n',pFold)

disp(foldTable)

fprintf('\nSuggested manuscript wording, if the improvement is meaningful:\n')
fprintf(['Ice thickness and temperature covary seasonally, so their global ' ...
    'correlation does not isolate the contribution of thickness. In ' ...
    'leave-one-dataset-out cross-validation, including thickness reduced ' ...
    'RMSE from %.1f to %.1f kg m-3 (Delta R2 = %.2f), indicating that ' ...
    'thickness captures density variability not explained by temperature ' ...
    'alone.\n'],cvTempCommon.rmse,cvFullCommon.rmse,deltaR2)

%% 6. Diagnostic figure
fig = figure('Units','inches','Position',[1 1 10.5 8.2]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% (a) Observed T-h relationship
ax1 = nexttile(tl,1);
scatter(ax1,T,h,38,month(date),'filled','MarkerEdgeColor','flat')
hold(ax1,'on')
ls = lsline(ax1);
set(ls,'LineWidth',1.2)
xlabel(ax1,'Core-average temperature T_i (°C)')
ylabel(ax1,'Sea-ice thickness h_i (m)')
title(ax1,sprintf('Observed association: r = %.2f',rPearson))
cb = colorbar(ax1);
cb.Ticks = 1:12;
cb.TickLabels = {'Jan','Feb','Mar','Apr','May','Jun', ...
                 'Jul','Aug','Sep','Oct','Nov','Dec'};
ylabel(cb,'Sampling month')
box(ax1,'on')
grid(ax1,'on')
text(ax1,-0.13,1.03,'(a)','Units','normalized','FontSize',10)

% (b) Cross-validated predictions
ax2 = nexttile(tl,2);
scatter(ax2,rho,cvTemp.yhat,36,'o','filled','MarkerFaceAlpha',0.65)
hold(ax2,'on')
scatter(ax2,rho,cvFull.yhat,36,'^','filled','MarkerFaceAlpha',0.65)
lims = [floor(min(rho)-2) ceil(max(rho)+2)];
plot(ax2,lims,lims,'k--','LineWidth',1)
xlim(ax2,lims)
ylim(ax2,lims)
axis(ax2,'square')
xlabel(ax2,'Observed density (kg m^{-3})')
ylabel(ax2,'Cross-validated density (kg m^{-3})')
legend(ax2, ...
    {sprintf('Temperature only, RMSE %.1f',cvTempCommon.rmse), ...
     sprintf('Temperature + thickness, RMSE %.1f',cvFullCommon.rmse), ...
     '1:1'}, ...
    'Location','northwest','Box','off')
box(ax2,'on')
grid(ax2,'on')
text(ax2,-0.13,1.03,'(b)','Units','normalized','FontSize',10)

% (c) Temperature-only residuals versus thickness
ax3 = nexttile(tl,3);
scatter(ax3,h(isCold),residTempCV(isCold),38,'o','filled','MarkerFaceAlpha',0.65)
hold(ax3,'on')
scatter(ax3,h(isWarm),residTempCV(isWarm),38,'^','filled','MarkerFaceAlpha',0.65)
yline(ax3,0,'k:','LineWidth',1)
add_regression_line(ax3,h(isCold),residTempCV(isCold),'--')
add_regression_line(ax3,h(isWarm),residTempCV(isWarm),'-')
xlabel(ax3,'Sea-ice thickness h_i (m)')
ylabel(ax3,'Temperature-only CV residual (kg m^{-3})')
legend(ax3, ...
    {sprintf('Cold: slope %.1f',coldStats.slope), ...
     sprintf('Warm: slope %.1f',warmStats.slope)}, ...
    'Location','best','Box','off')
box(ax3,'on')
grid(ax3,'on')
text(ax3,-0.13,1.03,'(c)','Units','normalized','FontSize',10)

% (d) Fold-level RMSE improvement
ax4 = nexttile(tl,4);
bar(ax4,categorical(foldTable.Dataset(validFolds)), ...
    foldTable.RMSE_temperature_only(validFolds) - ...
    foldTable.RMSE_full(validFolds))
yline(ax4,0,'k-','LineWidth',1)
ylabel(ax4,'RMSE improvement from thickness (kg m^{-3})')
xlabel(ax4,'Withheld dataset')
title(ax4,sprintf('Overall Delta RMSE = %.2f kg m^{-3}',deltaRMSE))
box(ax4,'on')
grid(ax4,'on')
xtickangle(ax4,35)
text(ax4,-0.13,1.03,'(d)','Units','normalized','FontSize',10)

set(findall(fig,'Type','axes'),'FontSize',10,'Toolbar',[])

if saveFigure
    outputFigure = fullfile(figureDir, ...
        'temperature_thickness_independence_diagnostic.png');
    exportgraphics(fig,outputFigure,'Resolution',300)
    fprintf('\nSaved figure to:\n%s\n',outputFigure)
end

%% Helpers
function fit = fit_temperature_only_model(rho,T,Tcrit_grid,min_n_cold,min_n_warm)
% Fit rho = wT*rhoCold + (1-wT)*rhoWarm without thickness.
% For each Tcrit, rhoCold and rhoWarm are estimated by least squares.

fit.rmse = inf;
fit.r2 = -inf;
sst = sum((rho-mean(rho,'omitnan')).^2,'omitnan');

for Tcrit = Tcrit_grid
    if sum(T <= Tcrit) < min_n_cold || sum(T > Tcrit) < min_n_warm
        continue
    end

    wT = temperature_weight(T,Tcrit);
    X = [wT, 1-wT];
    beta = X\rho;
    yhat = X*beta;
    rmse = sqrt(mean((rho-yhat).^2,'omitnan'));
    r2 = 1 - sum((rho-yhat).^2,'omitnan')/sst;

    if rmse < fit.rmse
        fit.rmse = rmse;
        fit.r2 = r2;
        fit.Tcrit = Tcrit;
        fit.rhoCold = beta(1);
        fit.rhoWarm = beta(2);
    end
end

if ~isfinite(fit.rmse)
    error('No valid temperature-only model found.')
end
end

function rhohat = predict_temperature_only_model(T,p)
Tcrit = p(1);
rhoCold = p(2);
rhoWarm = p(3);
wT = temperature_weight(T,Tcrit);
rhohat = wT.*rhoCold + (1-wT).*rhoWarm;
end

function cv = lodo_temperature_only_model(rho,T,dataset_id,p_full, ...
    Tcrit_grid,min_n_cold,min_n_warm)

groups = unique(dataset_id);
yhat = nan(size(rho));

for i = 1:numel(groups)
    idxTest = dataset_id == groups(i);
    idxTrain = ~idxTest;

    try
        fit_i = fit_temperature_only_model(rho(idxTrain),T(idxTrain), ...
            Tcrit_grid,min_n_cold,min_n_warm);
        p_i = [fit_i.Tcrit fit_i.rhoCold fit_i.rhoWarm];
    catch
        p_i = p_full;
    end

    yhat(idxTest) = predict_temperature_only_model(T(idxTest),p_i);
end

idx = isfinite(yhat);
cv = regression_metrics(rho(idx),yhat(idx),3);
cv.yhat = yhat;
end

function fit = fit_full_model(rho,T,h, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high, ...
    rho_cold_0_grid,Tcrit_grid,rho_warm_plateau_grid, ...
    hcrit_warm_grid,min_n_cold,min_n_warm)

fit.rmse = inf;
fit.r2 = -inf;
sst = sum((rho - mean(rho,'omitnan')).^2,'omitnan');

for iT = 1:numel(Tcrit_grid)
    Tcrit = Tcrit_grid(iT);

    if sum(T <= Tcrit) < min_n_cold || sum(T > Tcrit) < min_n_warm
        continue
    end

    wT = temperature_weight(T,Tcrit);

    for ih = 1:numel(hcrit_warm_grid)
        hcrit_warm = hcrit_warm_grid(ih);
        he_warm = min(h,hcrit_warm);
        Xh_warm = he_warm - hcrit_warm;

        for ic = 1:numel(rho_cold_0_grid)
            rho_cold_0 = rho_cold_0_grid(ic);
            rhoCold = cold_branch_density(h,rho_cold_0, ...
                h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);

            for ir = 1:numel(rho_warm_plateau_grid)
                rho_warm_plateau = rho_warm_plateau_grid(ir);
                base = wT.*rhoCold + (1-wT).*rho_warm_plateau;
                X = (1-wT).*Xh_warm;

                if all(abs(X) < eps)
                    continue
                end

                warm_slope = X\(rho-base);
                if ~isfinite(warm_slope) || warm_slope < 0
                    continue
                end

                yhat = base + X.*warm_slope;
                rmse = sqrt(mean((rho-yhat).^2,'omitnan'));
                r2 = 1 - sum((rho-yhat).^2,'omitnan')/sst;

                if rmse < fit.rmse
                    fit.rmse = rmse;
                    fit.r2 = r2;
                    fit.Tcrit = Tcrit;
                    fit.rho_warm_plateau = rho_warm_plateau;
                    fit.warm_slope = warm_slope;
                    fit.hcrit_warm = hcrit_warm;
                    fit.rho_cold_0 = rho_cold_0;
                end
            end
        end
    end
end

if ~isfinite(fit.rmse)
    error('No valid full model found.')
end
end

function rhohat = predict_full_model(T,h,p, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high)

Tcrit = p(1);
rho_warm_plateau = p(2);
warm_slope = p(3);
hcrit_warm = p(4);
rho_cold_0 = p(5);

rhoCold = cold_branch_density(h,rho_cold_0, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
he_warm = min(h,hcrit_warm);
rhoWarm = rho_warm_plateau + warm_slope.*(he_warm-hcrit_warm);
wT = temperature_weight(T,Tcrit);
rhohat = wT.*rhoCold + (1-wT).*rhoWarm;
end

function cv = lodo_full_model(rho,T,h,dataset_id,p_full, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high, ...
    rho_cold_0_grid,Tcrit_grid,rho_warm_plateau_grid, ...
    hcrit_warm_grid,min_n_cold,min_n_warm)

groups = unique(dataset_id);
yhat = nan(size(rho));

for i = 1:numel(groups)
    idxTest = dataset_id == groups(i);
    idxTrain = ~idxTest;

    try
        fit_i = fit_full_model(rho(idxTrain),T(idxTrain),h(idxTrain), ...
            h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high, ...
            rho_cold_0_grid,Tcrit_grid,rho_warm_plateau_grid, ...
            hcrit_warm_grid,min_n_cold,min_n_warm);
        p_i = [fit_i.Tcrit fit_i.rho_warm_plateau ...
               fit_i.warm_slope fit_i.hcrit_warm fit_i.rho_cold_0];
    catch
        p_i = p_full;
    end

    yhat(idxTest) = predict_full_model(T(idxTest),h(idxTest),p_i, ...
        h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
end

idx = isfinite(yhat);
cv = regression_metrics(rho(idx),yhat(idx),5);
cv.yhat = yhat;
end

function wT = temperature_weight(T,Tcrit)
wT = ones(size(T));
idxWarm = T > Tcrit;
wT(idxWarm) = T(idxWarm)./Tcrit;
wT = max(0,min(1,wT));
end

function rhoCold = cold_branch_density(h,rho_cold_0, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high)

rhoCold = nan(size(h));
idx1 = h <= h_cold_peak;
rhoCold(idx1) = rho_cold_0 + ...
    (rho_cold_peak-rho_cold_0).*h(idx1)./h_cold_peak;

idx2 = h > h_cold_peak;
rhoCold(idx2) = rho_cold_peak + ...
    (rho_cold_high-rho_cold_peak).* ...
    (min(h(idx2),h_cold_high)-h_cold_peak)./ ...
    (h_cold_high-h_cold_peak);
end

function metrics = regression_metrics(y,yhat,k)
idx = isfinite(y) & isfinite(yhat);
y = y(idx);
yhat = yhat(idx);
e = y-yhat;
n = numel(y);
sse = sum(e.^2,'omitnan');
sst = sum((y-mean(y,'omitnan')).^2,'omitnan');
metrics.rmse = sqrt(mean(e.^2,'omitnan'));
metrics.r2 = 1-sse/sst;
metrics.r2adj = 1-(1-metrics.r2)*(n-1)/max(n-k-1,1);
metrics.aic = n*log(sse/n)+2*k;
metrics.bic = n*log(sse/n)+k*log(n);
metrics.N = n;
end

function stats = residual_thickness_test(h,resid)
idx = isfinite(h) & isfinite(resid);
h = h(idx);
resid = resid(idx);
stats.N = numel(h);

if stats.N < 3 || range(h) == 0
    stats.slope = NaN;
    stats.intercept = NaN;
    stats.r = NaN;
    stats.p = NaN;
    return
end

b = [ones(size(h)) h]\resid;
stats.intercept = b(1);
stats.slope = b(2);
[stats.r,stats.p] = corr(h,resid,'Type','Pearson','Rows','complete');
end

function foldTable = build_fold_table(rho,dataset_id,yhatTemp,yhatFull)
groups = unique(dataset_id);
Dataset = strings(numel(groups),1);
N = zeros(numel(groups),1);
RMSE_temperature_only = nan(numel(groups),1);
RMSE_full = nan(numel(groups),1);

for i = 1:numel(groups)
    idx = dataset_id == groups(i) & ...
          isfinite(yhatTemp) & isfinite(yhatFull);
    Dataset(i) = groups(i);
    N(i) = sum(idx);
    if N(i) > 0
        RMSE_temperature_only(i) = sqrt(mean( ...
            (rho(idx)-yhatTemp(idx)).^2,'omitnan'));
        RMSE_full(i) = sqrt(mean( ...
            (rho(idx)-yhatFull(idx)).^2,'omitnan'));
    end
end

Delta_RMSE = RMSE_temperature_only-RMSE_full;
foldTable = table(Dataset,N,RMSE_temperature_only,RMSE_full,Delta_RMSE);
end

function add_regression_line(ax,x,y,lineStyle)
idx = isfinite(x) & isfinite(y);
if sum(idx) < 3 || range(x(idx)) == 0
    return
end
b = [ones(sum(idx),1) x(idx)]\y(idx);
xlineData = linspace(min(x(idx)),max(x(idx)),100)';
ylineData = b(1)+b(2)*xlineData;
plot(ax,xlineData,ylineData,lineStyle,'LineWidth',1.4, ...
    'HandleVisibility','off')
end
