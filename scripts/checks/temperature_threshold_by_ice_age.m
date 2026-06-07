clear
close all
clc

repoRoot = 'C:\Users\evsalg001\Documents\MATLAB\Density parametrization\script';

finalDir = fullfile(repoRoot, ...
    'data', 'final');

figureDir = fullfile(repoRoot, ...
    'figures');

colormapDir = fullfile(repoRoot, ...
    'data', 'colormaps');

if ~exist(figureDir, 'dir')
    mkdir(figureDir)
end

inputFile = fullfile(finalDir, ...
    'snow_model_matchups_with_models.mat');

load(inputFile,'D')

S = load(fullfile(colormapDir,'lipari.mat'));
lipari = S.lipari;

densityMode = 1; % 1 = lab, 2 = connected, 3 = disconnected
rho_cols = {"density_lab_kgm3", ...
            "density_insitu_connected_kgm3", ...
            "density_insitu_disconnected_kgm3"};
rho_col = rho_cols{densityMode};

requiredVars = ["date","ice_age","dataset","temperature_C","ice_thickness_m",rho_col];
missingVars = requiredVars(~ismember(requiredVars,string(D.Properties.VariableNames)));
if ~isempty(missingVars)
    error('D is missing required variables: %s',strjoin(cellstr(missingVars),', '))
end

% Cold branch for T_i <= Tcrit:
% bilinear in ice thickness with an independent peak near 1 m.
h_cold_peak = 1.10;
rho_cold_peak = 911;
h_cold_high = 3.50;
rho_cold_high = 908.5;
rho_cold_0_grid = 905:1:912;

% Warm branch for T_i > Tcrit:
% bilinear in ice thickness with its own critical thickness near 2 m.
Tcrit_grid = -3.5:0.05:-1.5;
hcrit_warm_grid = 1.8:0.1:2.3;
rho_warm_plateau_grid = 895:1:905;

min_n_cold = 8;
min_n_warm = 8;

% Bootstrap is used only for the confidence bands in the two-panel figure.
doBootstrap = true;
nboot = 200;
rng(1)

if ~isdatetime(D.date)
    D.date = datetime(D.date);
end

if ismember("is_freeze",string(D.Properties.VariableNames))
    is_freeze_all = logical(D.is_freeze);
else
    is_freeze_all = false(height(D),1);
end

t = D.date;
t.Year = 2020;
T = D.temperature_C;
rho = D.(rho_col);
h = D.ice_thickness_m;
dataset_id = string(D.dataset);
ice_type = upper(string(D.ice_age));
is_freeze = is_freeze_all;
is_svalbard = contains(lower(dataset_id),"svalbard");

valid = ~isnat(t) & isfinite(T) & isfinite(rho) & isfinite(h) & T < 0 & h > 0;

t = t(valid);
T = T(valid);
rho = rho(valid);
h = h(valid);
dataset_id = dataset_id(valid);
ice_type = ice_type(valid);
is_freeze = is_freeze(valid);
is_svalbard = is_svalbard(valid);

is_fyi_all = ice_type == "FYI" & ~is_svalbard & ~is_freeze;
is_old_all = ice_type ~= "FYI" & ~is_svalbard & ~is_freeze;
is_sva_all = is_svalbard & ~is_freeze;

T_fit = T(~is_freeze);
rho_fit = rho(~is_freeze);
h_fit = h(~is_freeze);
dataset_fit = dataset_id(~is_freeze);
N = numel(rho_fit);

fit = fit_fast_fixed_hcrit_model( ...
    rho_fit,T_fit,h_fit, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid, ...
    Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm);

p = [fit.Tcrit fit.rho_warm_plateau fit.warm_slope fit.hcrit_warm fit.rho_cold_0];

rhohat = predict_fast_fixed_hcrit_model( ...
    T_fit,h_fit,p,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);

metrics = regression_metrics(rho_fit,rhohat,numel(p));

if doBootstrap
    B_fit = bootstrap_fast_fixed_hcrit_model( ...
        rho_fit,T_fit,h_fit,dataset_fit,nboot, ...
        h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid, ...
        Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm);
else
    B_fit = repmat(p,nboot,1);
end

h_sel = [0.5 1.5 2.0];
Tplot = linspace(-10,-0.05,400)';
rho_sel_1 = predict_fast_fixed_hcrit_model(Tplot,h_sel(1)*ones(size(Tplot)),p,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
rho_sel_2 = predict_fast_fixed_hcrit_model(Tplot,h_sel(2)*ones(size(Tplot)),p,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
rho_sel_3 = predict_fast_fixed_hcrit_model(Tplot,h_sel(3)*ones(size(Tplot)),p,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);

T_sel = [p(1) -1.5 -0.5];
hplot = linspace(0,3.5,400)';
rho_T1 = predict_fast_fixed_hcrit_model(T_sel(1)*ones(size(hplot)),hplot,p,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
rho_T2 = predict_fast_fixed_hcrit_model(T_sel(2)*ones(size(hplot)),hplot,p,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
rho_T3 = predict_fast_fixed_hcrit_model(T_sel(3)*ones(size(hplot)),hplot,p,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);

fs_ax = 10.5;
fs_txt = 9.0;
lw = 1.2;
ms = 30;
ms_ca_fyi = 6.0;
ms_ca_syi = 6.0;
ms_sva = 6.2;
lw_open = 1.0;
Tmin_time = -8;
Tmax_time = 0;

fig1 = figure;
set(fig1,'Units','inches','Position',[1.1 1.0 10.5 4.4]);
tl = tiledlayout(fig1,1,2,'TileSpacing','compact','Padding','compact');

ax2 = nexttile(tl,1);
hold(ax2,'on')
scatter(ax2,T(is_fyi_all),rho(is_fyi_all),ms,h(is_fyi_all),'filled','MarkerEdgeColor','flat','LineWidth',0.35)
scatter(ax2,T(is_old_all),rho(is_old_all),ms,h(is_old_all),'MarkerFaceColor','none','MarkerEdgeColor','flat','LineWidth',1.0)
scatter(ax2,T(is_sva_all),rho(is_sva_all),ms,h(is_sva_all),'>','filled','MarkerEdgeColor','flat','LineWidth',0.35)
colormap(ax2,lipari)
clim(ax2,[min(h,[],'omitnan') max(h,[],'omitnan')])
cb2 = colorbar(ax2);
ylabel(cb2,'Sea-ice thickness (m)','FontSize',fs_ax)
cmap = colormap(ax2);
cl = ax2.CLim;
ncol = size(cmap,1);
get_color = @(val) cmap(max(1,min(ncol,round((val-cl(1))/(cl(2)-cl(1))*(ncol-1))+1)),:);
c1p = get_color(h_sel(1));
c2p = get_color(h_sel(2));
c3p = get_color(h_sel(3));

h_ci2 = add_bootstrap_ci_band(ax2,Tplot,h_sel(1)*ones(size(Tplot)),B_fit,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
set(h_ci2,'DisplayName','95% bootstrap CI','HandleVisibility','on')
add_bootstrap_ci_band(ax2,Tplot,h_sel(2)*ones(size(Tplot)),B_fit,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
add_bootstrap_ci_band(ax2,Tplot,h_sel(3)*ones(size(Tplot)),B_fit,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);

p_h1 = plot(ax2,Tplot,rho_sel_1,'-','Color',c1p,'LineWidth',lw);
p_h2 = plot(ax2,Tplot,rho_sel_2,'-','Color',c2p,'LineWidth',lw);
p_h3 = plot(ax2,Tplot,rho_sel_3,'-','Color',c3p,'LineWidth',lw);
p_fyi2 = plot(ax2,nan,nan,'o','MarkerSize',ms_ca_fyi,'MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','none');
p_syi2 = plot(ax2,nan,nan,'o','MarkerSize',ms_ca_syi,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',lw_open,'LineStyle','none');
p_sva2 = plot(ax2,nan,nan,'>','MarkerSize',ms_sva,'MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','none');
leg2 = legend(ax2,[p_fyi2 p_syi2 p_sva2 p_h1 p_h2 p_h3 h_ci2], ...
    {'FYI','SYI \& MYI','Svalbard', ...
    sprintf('$h_i = %.1f\\ \\mathrm{m}$',h_sel(1)), ...
    sprintf('$h_i = %.1f\\ \\mathrm{m}$',h_sel(2)), ...
    sprintf('$h_i = %.1f\\ \\mathrm{m}$',h_sel(3)), ...
    '95\% bootstrap CI'}, ...
    'Location','west', ...
    'Box','off', ...
    'FontSize',fs_txt, ...
    'Interpreter','latex');
leg2.ItemTokenSize = [14 4];

eqnTxt = sprintf([ ...
    '$w_T=T_i/(%.1f),\\quad T_i>%.1f^\\circ\\mathrm{C}$' newline ...
    '$w_T=1,\\quad T_i\\leq%.1f^\\circ\\mathrm{C}$' newline ...
    '$\\rho_w=%.0f%+.0f\\min(h_i,%.1f\\ \\mathrm{m})$' newline ...
    '$\\rho_c=%.0f+(%.0f-%.0f)h_i/%.1f,\\quad h_i\\leq%.1f\\ \\mathrm{m}$' newline ...
    '$\\rho_c=%.0f%+.1f(h_i-%.1f),\\quad h_i>%.1f\\ \\mathrm{m}$' newline ...
    '$\\rho=w_T\\rho_c+(1-w_T)\\rho_w$' newline ...
    '$R^2=%.2f,\\ RMSE=%.1f\\ \\mathrm{kg\\ m^{-3}},\\ N=%d$'], ...
    p(1),p(1), ...
    p(1), ...
    round(p(2)-p(3)*p(4)), p(3), p(4), ...
    p(5),rho_cold_peak,p(5),h_cold_peak,h_cold_peak, ...
    rho_cold_peak,(rho_cold_high-rho_cold_peak)/(h_cold_high-h_cold_peak),h_cold_peak,h_cold_peak, ...
    metrics.r2,metrics.rmse,N);
text(ax2,0.02,0.01,eqnTxt,'Units','normalized','Interpreter','latex','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',fs_txt + 0.8)
xlabel(ax2,'Sea-ice temperature T_i (°C)','Interpreter','tex')
ylabel(ax2,'Sea-ice density \rho (kg m^{-3})','Interpreter','tex')
set(ax2,'FontSize',fs_ax,'FontWeight','normal')
xlim(ax2,[-10 0])
ylim(ax2,[850 920])
box(ax2,'on')
text(ax2,-0.145,1.02,'(a)','Units','normalized','FontSize',10,'FontWeight','normal')

ax3 = nexttile(tl,2);
hold(ax3,'on')
scatter(ax3,h(is_fyi_all),rho(is_fyi_all),ms,T(is_fyi_all),'filled','MarkerEdgeColor','flat','LineWidth',0.35)
scatter(ax3,h(is_old_all),rho(is_old_all),ms,T(is_old_all),'MarkerFaceColor','none','MarkerEdgeColor','flat','LineWidth',1.0)
scatter(ax3,h(is_sva_all),rho(is_sva_all),ms,T(is_sva_all),'>','filled','MarkerEdgeColor','flat','LineWidth',0.35)
colormap(ax3,lipari)
clim(ax3,[Tmin_time Tmax_time])
cb3 = colorbar(ax3);
ylabel(cb3,'Sea-ice temperature (°C)','FontSize',fs_ax)
cmap3 = colormap(ax3);
cl3 = ax3.CLim;
ncol3 = size(cmap3,1);
get_color_T = @(val) cmap3(max(1,min(ncol3,round((val-cl3(1))/(cl3(2)-cl3(1))*(ncol3-1))+1)),:);
cT1 = get_color_T(T_sel(1));
cT2 = get_color_T(T_sel(2));
cT3 = get_color_T(T_sel(3));

h_ci3 = add_bootstrap_ci_band(ax3,T_sel(1)*ones(size(hplot)),hplot,B_fit,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
set(h_ci3,'DisplayName','95% bootstrap CI','HandleVisibility','on')
add_bootstrap_ci_band(ax3,T_sel(2)*ones(size(hplot)),hplot,B_fit,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
add_bootstrap_ci_band(ax3,T_sel(3)*ones(size(hplot)),hplot,B_fit,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);

p_T1 = plot(ax3,hplot,rho_T1,'-','Color',cT1,'LineWidth',lw);
p_T2 = plot(ax3,hplot,rho_T2,'-','Color',cT2,'LineWidth',lw);
p_T3 = plot(ax3,hplot,rho_T3,'-','Color',cT3,'LineWidth',lw);

p_fyi3 = plot(ax3,nan,nan,'o','MarkerSize',ms_ca_fyi,'MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','none');
p_syi3 = plot(ax3,nan,nan,'o','MarkerSize',ms_ca_syi,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',lw_open,'LineStyle','none');
p_sva3 = plot(ax3,nan,nan,'>','MarkerSize',ms_sva,'MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','none');

leg3 = legend(ax3, ...
    [p_fyi3 p_syi3 p_sva3 p_T1 p_T2 p_T3 h_ci3], ...
    {'FYI','SYI \& MYI','Svalbard', ...
     sprintf('$T_i = %.1f^{\\circ}\\mathrm{C}$',T_sel(1)), ...
     sprintf('$T_i = %.1f^{\\circ}\\mathrm{C}$',T_sel(2)), ...
     sprintf('$T_i = %.1f^{\\circ}\\mathrm{C}$',T_sel(3)), ...
     '95\% bootstrap CI'}, ...
    'Location','southeast', ...
    'Box','off', ...
    'FontSize',fs_txt, ...
    'Interpreter','latex');
leg3.AutoUpdate = 'off';
leg3.ItemTokenSize = [14 4];
xlabel(ax3,'Sea-ice thickness h_i (m)','Interpreter','tex')
ylabel(ax3,'Sea-ice density \rho (kg m^{-3})','Interpreter','tex')
set(ax3,'FontSize',fs_ax,'FontWeight','normal')
xlim(ax3,[0 3.5])
ylim(ax3,[850 920])
box(ax3,'on')
text(ax3,-0.145,1.02,'(b)','Units','normalized','FontSize',10,'FontWeight','normal')

% outputFigure = fullfile(figureDir,'Fig_temperature_density_two_panel.png');
% set(findall(fig1,'Type','axes'),'Toolbar',[])
% exportgraphics(fig1,outputFigure,'Resolution',300)
% fprintf('Generated two-panel temperature-density figure using %s density (%d observations).\n', rho_col, N)
% fprintf('Saved figure to:\n%s\n', outputFigure)

%% Ice-type threshold sensitivity test (FYI vs SYI/MYI)
% This test addresses whether the temperature transition can be constrained
% separately for FYI and old ice (SYI/MYI). It should be interpreted as a
% sensitivity analysis, not as a replacement for the pooled parameterization.

fprintf('\nIce-type threshold sensitivity test:\n')

is_fyi_fit = ice_type(~is_freeze) == "FYI" & ~is_svalbard(~is_freeze);
is_old_fit = ice_type(~is_freeze) ~= "FYI" & ~is_svalbard(~is_freeze);

fprintf('FYI non-Svalbard N = %d; SYI/MYI non-Svalbard N = %d\n',sum(is_fyi_fit),sum(is_old_fit))

% 1) Fit the full model separately to FYI and SYI/MYI where possible.
fit_fyi = [];
fit_old = [];

if sum(is_fyi_fit) >= 25
    fit_fyi = fit_fast_fixed_hcrit_model( ...
        rho_fit(is_fyi_fit),T_fit(is_fyi_fit),h_fit(is_fyi_fit), ...
        h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid, ...
        Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm);
    p_fyi = [fit_fyi.Tcrit fit_fyi.rho_warm_plateau fit_fyi.warm_slope fit_fyi.hcrit_warm fit_fyi.rho_cold_0];
    yhat_fyi = predict_fast_fixed_hcrit_model(T_fit(is_fyi_fit),h_fit(is_fyi_fit),p_fyi, ...
        h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
    met_fyi = regression_metrics(rho_fit(is_fyi_fit),yhat_fyi,numel(p_fyi));
    fprintf('FYI-only full refit: Tcrit = %.2f C, R2 = %.2f, RMSE = %.1f kg m-3\n', ...
        p_fyi(1),met_fyi.r2,met_fyi.rmse)
else
    fprintf('FYI-only full refit skipped: too few observations.\n')
end

if sum(is_old_fit) >= 25 && sum(T_fit(is_old_fit) <= min(Tcrit_grid)) >= min_n_cold && sum(T_fit(is_old_fit) > max(Tcrit_grid)) >= min_n_warm
    fit_old = fit_fast_fixed_hcrit_model( ...
        rho_fit(is_old_fit),T_fit(is_old_fit),h_fit(is_old_fit), ...
        h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid, ...
        Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm);
    p_old = [fit_old.Tcrit fit_old.rho_warm_plateau fit_old.warm_slope fit_old.hcrit_warm fit_old.rho_cold_0];
    yhat_old = predict_fast_fixed_hcrit_model(T_fit(is_old_fit),h_fit(is_old_fit),p_old, ...
        h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
    met_old = regression_metrics(rho_fit(is_old_fit),yhat_old,numel(p_old));
    fprintf('SYI/MYI-only full refit: Tcrit = %.2f C, R2 = %.2f, RMSE = %.1f kg m-3\n', ...
        p_old(1),met_old.r2,met_old.rmse)
else
    fprintf('SYI/MYI-only full refit skipped or weakly constrained: too few cold/warm observations for a stable 5-parameter fit.\n')
end

% 2) More conservative test: keep all non-threshold parameters fixed from
% the pooled model and scan only Tcrit separately for FYI and old ice.
% This asks whether the data support a meaningful Tcrit split without giving
% each ice type a full independent 5-parameter model.
non_svalbard_fit = ~is_svalbard(~is_freeze);
T_test = T_fit(non_svalbard_fit);
h_test = h_fit(non_svalbard_fit);
rho_test = rho_fit(non_svalbard_fit);
ice_test = ice_type(~is_freeze);
ice_test = ice_test(non_svalbard_fit);
is_fyi_test = ice_test == "FYI";
is_old_test = ice_test ~= "FYI";

p_base = p;

% Common-threshold SSE using the pooled model evaluated on non-Svalbard FYI/SYI/MYI cores.
yhat_common = predict_fast_fixed_hcrit_model(T_test,h_test,p_base, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
sse_common = sum((rho_test - yhat_common).^2,'omitnan');
rmse_common = sqrt(mean((rho_test - yhat_common).^2,'omitnan'));

sse_split = nan(numel(Tcrit_grid),numel(Tcrit_grid));

for i = 1:numel(Tcrit_grid)
    for j = 1:numel(Tcrit_grid)
        p_tmp_fyi = p_base;
        p_tmp_old = p_base;
        p_tmp_fyi(1) = Tcrit_grid(i);
        p_tmp_old(1) = Tcrit_grid(j);

        yhat_tmp = nan(size(rho_test));
        yhat_tmp(is_fyi_test) = predict_fast_fixed_hcrit_model(T_test(is_fyi_test),h_test(is_fyi_test),p_tmp_fyi, ...
            h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
        yhat_tmp(is_old_test) = predict_fast_fixed_hcrit_model(T_test(is_old_test),h_test(is_old_test),p_tmp_old, ...
            h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);

        % Require both sides of each threshold to be represented when possible.
        if sum(T_test(is_fyi_test) <= Tcrit_grid(i)) < min_n_cold || sum(T_test(is_fyi_test) > Tcrit_grid(i)) < min_n_warm
            continue
        end
        if sum(T_test(is_old_test) <= Tcrit_grid(j)) < 3 || sum(T_test(is_old_test) > Tcrit_grid(j)) < 3
            continue
        end

        sse_split(i,j) = sum((rho_test - yhat_tmp).^2,'omitnan');
    end
end

[min_sse_split,idx_min] = min(sse_split(:),[],'omitnan');
[i_best,j_best] = ind2sub(size(sse_split),idx_min);
Tcrit_fyi_best = Tcrit_grid(i_best);
Tcrit_old_best = Tcrit_grid(j_best);
rmse_split = sqrt(min_sse_split/numel(rho_test));

% Compare by AIC/BIC. Split-threshold model has one extra parameter.
n_test = numel(rho_test);
k_common = numel(p_base);
k_split = numel(p_base) + 1;
aic_common = n_test*log(sse_common/n_test) + 2*k_common;
aic_split = n_test*log(min_sse_split/n_test) + 2*k_split;
bic_common = n_test*log(sse_common/n_test) + k_common*log(n_test);
bic_split = n_test*log(min_sse_split/n_test) + k_split*log(n_test);

fprintf('Common Tcrit model: Tcrit = %.2f C, RMSE = %.2f kg m-3, AIC = %.1f, BIC = %.1f\n', ...
    p_base(1),rmse_common,aic_common,bic_common)
fprintf('Split Tcrit scan: FYI Tcrit = %.2f C, SYI/MYI Tcrit = %.2f C, RMSE = %.2f kg m-3, AIC = %.1f, BIC = %.1f\n', ...
    Tcrit_fyi_best,Tcrit_old_best,rmse_split,aic_split,bic_split)
fprintf('Delta RMSE common - split = %.2f kg m-3; Delta AIC split-common = %.1f; Delta BIC split-common = %.1f\n', ...
    rmse_common-rmse_split,aic_split-aic_common,bic_split-bic_common)

% 3) Profile curves for Tcrit by ice type with other parameters fixed.
rmse_fyi_profile = nan(size(Tcrit_grid));
rmse_old_profile = nan(size(Tcrit_grid));

for i = 1:numel(Tcrit_grid)
    p_tmp = p_base;
    p_tmp(1) = Tcrit_grid(i);
    if sum(T_test(is_fyi_test) <= Tcrit_grid(i)) >= min_n_cold && sum(T_test(is_fyi_test) > Tcrit_grid(i)) >= min_n_warm
        yhat_tmp = predict_fast_fixed_hcrit_model(T_test(is_fyi_test),h_test(is_fyi_test),p_tmp, ...
            h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
        rmse_fyi_profile(i) = sqrt(mean((rho_test(is_fyi_test) - yhat_tmp).^2,'omitnan'));
    end
    if sum(T_test(is_old_test) <= Tcrit_grid(i)) >= 3 && sum(T_test(is_old_test) > Tcrit_grid(i)) >= 3
        yhat_tmp = predict_fast_fixed_hcrit_model(T_test(is_old_test),h_test(is_old_test),p_tmp, ...
            h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
        rmse_old_profile(i) = sqrt(mean((rho_test(is_old_test) - yhat_tmp).^2,'omitnan'));
    end
end

fig_thr = figure;
set(fig_thr,'Units','inches','Position',[1.2 1.0 6.0 4.0]);
hold on
plot(Tcrit_grid,rmse_fyi_profile,'-','LineWidth',1.4)
plot(Tcrit_grid,rmse_old_profile,'--','LineWidth',1.4)
xline(p_base(1),':','LineWidth',1.0)
xlabel('Assumed transition temperature T_{crit} (°C)')
ylabel('RMSE (kg m^{-3})')
legend({'FYI','SYI/MYI','pooled T_{crit}'},'Location','best','Box','off')
box on
grid on
% outputThresholdFigure = fullfile(figureDir,'Tcrit_profile_by_ice_type.png');
% exportgraphics(fig_thr,outputThresholdFigure,'Resolution',300)
% fprintf('Saved threshold profile figure to:\n%s\n',outputThresholdFigure)

%% Helpers
function fit = fit_fast_fixed_hcrit_model(rho,T,h,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid,Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm)

fit.rmse = inf;
fit.r2 = -inf;
sst = sum((rho - mean(rho,'omitnan')).^2,'omitnan');

for iT = 1:numel(Tcrit_grid)
    Tcrit = Tcrit_grid(iT);

    if sum(T <= Tcrit) < min_n_cold || sum(T > Tcrit) < min_n_warm
        continue
    end

    wT = ones(size(T));
    idxWarm = T > Tcrit;
    wT(idxWarm) = T(idxWarm)./Tcrit;
    wT = max(0,min(1,wT));

    for ih = 1:numel(hcrit_warm_grid)
        hcrit_warm = hcrit_warm_grid(ih);
        he_warm = min(h,hcrit_warm);
        Xh_warm = he_warm - hcrit_warm;

        for ic = 1:numel(rho_cold_0_grid)
            rho_cold_0 = rho_cold_0_grid(ic);
            rhoCold = cold_branch_density(h,rho_cold_0,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);

            for ir = 1:numel(rho_warm_plateau_grid)
                rho_warm_plateau = rho_warm_plateau_grid(ir);
                base = wT.*rhoCold + (1 - wT).*rho_warm_plateau;
                X = (1 - wT).*Xh_warm;

                if all(abs(X) < eps)
                    continue
                end

                warm_slope = X\(rho - base);

                if ~isfinite(warm_slope) || warm_slope < 0
                    continue
                end

                yhat = base + X.*warm_slope;
                rmse = sqrt(mean((rho - yhat).^2,'omitnan'));
                r2 = 1 - sum((rho - yhat).^2,'omitnan')/sst;

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
    error('No valid model found. Widen Tcrit_grid, hcrit_warm_grid, rho_cold_0_grid, or rho_warm_plateau_grid.')
end
end

function rhohat = predict_fast_fixed_hcrit_model(T,h,p,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high)

Tcrit = p(1);
rho_warm_plateau = p(2);
warm_slope = p(3);
hcrit_warm = p(4);
rho_cold_0 = p(5);

rhoCold = cold_branch_density(h,rho_cold_0,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
he_warm = min(h,hcrit_warm);
rhoWarm = rho_warm_plateau + warm_slope.*(he_warm - hcrit_warm);

wT = ones(size(T));
idxWarm = T > Tcrit;
wT(idxWarm) = T(idxWarm)./Tcrit;
wT = max(0,min(1,wT));

rhohat = wT.*rhoCold + (1 - wT).*rhoWarm;
end

function rhoCold = cold_branch_density(h,rho_cold_0,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high)

% Bilinear cold branch, independent from the warm-branch critical thickness:
% segment 1: increasing density from h_i = 0 to h_cold_peak (~1 m)
% segment 2: slowly decreasing density for h_i > h_cold_peak
rhoCold = nan(size(h));

idx1 = h <= h_cold_peak;
rhoCold(idx1) = rho_cold_0 + ...
    (rho_cold_peak - rho_cold_0).*h(idx1)./h_cold_peak;

idx2 = h > h_cold_peak;
rhoCold(idx2) = rho_cold_peak + ...
    (rho_cold_high - rho_cold_peak).* ...
    (min(h(idx2),h_cold_high) - h_cold_peak)./(h_cold_high - h_cold_peak);
end

function cv = leave_one_dataset_out_cv_fast_fixed_hcrit(rho,T,h,dataset_id,p_full,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid,Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm)

groups = unique(dataset_id);
yhat_all = nan(size(rho));

for i = 1:numel(groups)
    idxTest = dataset_id == groups(i);
    idxTrain = ~idxTest;

    if sum(idxTrain) < 25 || sum(idxTest) < 1
        continue
    end

    try
        fit_i = fit_fast_fixed_hcrit_model(rho(idxTrain),T(idxTrain),h(idxTrain),h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid,Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm);
        p_i = [fit_i.Tcrit fit_i.rho_warm_plateau fit_i.warm_slope fit_i.hcrit_warm fit_i.rho_cold_0];
    catch
        p_i = p_full;
    end

    yhat_all(idxTest) = predict_fast_fixed_hcrit_model(T(idxTest),h(idxTest),p_i,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
end

idx = isfinite(yhat_all);
cv = regression_metrics(rho(idx),yhat_all(idx),5);
cv.yhat = yhat_all;
end

function metrics = regression_metrics(y,yhat,k)

e = y - yhat;
n = numel(y);
sse = sum(e.^2,'omitnan');
sst = sum((y - mean(y,'omitnan')).^2,'omitnan');
metrics.rmse = sqrt(mean(e.^2,'omitnan'));
metrics.r2 = 1 - sse/sst;
metrics.r2adj = 1 - (1 - metrics.r2)*(n - 1)/max(n - k - 1,1);
metrics.aic = n*log(sse/n) + 2*k;
metrics.bic = n*log(sse/n) + k*log(n);
end

function B = bootstrap_fast_fixed_hcrit_model(rho,T,h,dataset_id,nboot,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid,Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm)

n = numel(rho);
B = nan(nboot,5);
groups = unique(dataset_id);
G = numel(groups);

for b = 1:nboot
    pickGroups = groups(randi(G,G,1));
    idx = false(n,1);

    for j = 1:numel(pickGroups)
        idx = idx | dataset_id == pickGroups(j);
    end

    if sum(idx) < 25
        idx = randi(n,n,1);
    end

    try
        fit_b = fit_fast_fixed_hcrit_model(rho(idx),T(idx),h(idx),h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid,Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm);
        B(b,:) = [fit_b.Tcrit fit_b.rho_warm_plateau fit_b.warm_slope fit_b.hcrit_warm fit_b.rho_cold_0];
    catch
        B(b,:) = nan(1,5);
    end
end

B = B(all(isfinite(B),2),:);
end

function hband = add_bootstrap_ci_band(ax,T_grid,h_grid,B,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high)

rho_boot = nan(numel(T_grid),size(B,1));

for b = 1:size(B,1)
    rho_boot(:,b) = predict_fast_fixed_hcrit_model(T_grid,h_grid,B(b,:),h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);
end

rho_lo = prctile(rho_boot,2.5,2);
rho_hi = prctile(rho_boot,97.5,2);
x = T_grid;

if range(T_grid) == 0
    x = h_grid;
end

hband = fill(ax,[x(:); flipud(x(:))],[rho_lo(:); flipud(rho_hi(:))], ...
    [0.75 0.75 0.75],'FaceAlpha',0.20,'EdgeColor','none','HandleVisibility','off');
end
