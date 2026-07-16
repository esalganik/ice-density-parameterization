function analyze_density_model(repoRoot)
%
% ANALYZE_DENSITY_MODEL
%
% Derives the observation-based Arctic sea-ice density parameterization
% presented in Figure 1 of Salganik et al.
%
% Workflow:
%   1. Load merged observational density and SnowModel matchups.
%   2. Select density definition (laboratory or in-situ).
%   3. Fit a temperature-thickness density model.
%   4. Estimate parameter uncertainty using bootstrap resampling.
%   5. Evaluate model skill using leave-one-dataset-out cross-validation.
%   6. Generate Figure 1.
%
% Input:
%   data/final/snow_model_matchups_with_models.mat
%
% Output:
%   figures/Fig1.png
%
% Notes:
%   The fitted model combines a cold-regime thickness-dependent density
%   parameterization with a warm-regime density reduction controlled by
%   sea-ice temperature.
%

close all

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

densityMode = 1;
% Density definition:
% 1 = laboratory measurements
% 2 = in-situ connected porosity
% 3 = in-situ disconnected porosity
rho_cols = {"density_lab_kgm3", ...
            "density_insitu_connected_kgm3", ...
            "density_insitu_disconnected_kgm3"};
rho_names = {"lab","connected","disconnected"};
rho_col = rho_cols{densityMode};

requiredVars = ["date","ice_age","dataset","temperature_C","ice_thickness_m",rho_col];
missingVars = requiredVars(~ismember(requiredVars,string(D.Properties.VariableNames)));
if ~isempty(missingVars)
    error('D is missing required variables: %s',strjoin(cellstr(missingVars),', '))
end

% Cold-regime density parameterization (T_i <= T_crit).
% Density varies piecewise linearly with ice thickness, reaching a maximum
% at h_cold_peak before decreasing toward thicker ice.
h_cold_peak = 1.10;
rho_cold_peak = 911;
h_cold_high = 3.50;
rho_cold_high = 908.5;
rho_cold_0_grid = 905:1:912;

% Warm-regime density parameterization (T_i > T_crit).
% Density varies linearly with effective thickness up to hcrit_warm and approaches a plateau for thicker ice.
% Candidate parameter values explored during grid-search optimization.
Tcrit_grid = -3.0:0.05:-2.0;
hcrit_warm_grid = 1.8:0.1:2.3;
rho_warm_plateau_grid = 895:1:905;

min_n_cold = 8;
min_n_warm = 8;

doBootstrap = true;
doLODO = true; % Optional diagnostic: leave-one-dataset-out cross-validation.
nboot = 1000;
rng(1)

if ~isdatetime(D.date)
    D.date = datetime(D.date);
end

t = D.date;
% Replace calendar year with a common reference year so observations from multiple campaigns can be combined into a climatological seasonal cycle.
t.Year = 2020;
T = D.temperature_C;
rho = D.(rho_col);
h = D.ice_thickness_m;
dataset_id = string(D.dataset);
ice_type = upper(string(D.ice_age));
is_svalbard = contains(lower(dataset_id),"svalbard");

valid = ~isnat(t) & isfinite(T) & isfinite(rho) & isfinite(h) & T < 0 & h > 0;

t = t(valid);
T = T(valid);
rho = rho(valid);
h = h(valid);
dataset_id = dataset_id(valid);
ice_type = ice_type(valid);
is_svalbard = is_svalbard(valid);

% Separate Arctic FYI, older ice (SYI/MYI), and Svalbard fjord ice for
% visualization. Svalbard observations are shown separately because they
% represent a coastal/fjord environment rather than drifting pack ice.
is_fyi_all = ice_type == "FYI" & ~is_svalbard;
is_old_all = ice_type ~= "FYI" & ~is_svalbard;
is_sva_all = is_svalbard;

T_fit = T;
rho_fit = rho;
h_fit = h;
dataset_fit = dataset_id;
N = numel(rho_fit);

% Fit a hybrid density model:
%   - Cold regime: piecewise-linear dependence on ice thickness.
%   - Warm regime: temperature-dependent transition toward a lower-density
%     thickness-controlled branch.

fit = fit_fast_fixed_hcrit_model( ...
    rho_fit,T_fit,h_fit, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid, ...
    Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm);

p = [fit.Tcrit fit.rho_warm_plateau fit.warm_slope fit.hcrit_warm fit.rho_cold_0];

rhohat = predict_fast_fixed_hcrit_model( ...
    T_fit,h_fit,p,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);

metrics = regression_metrics(rho_fit,rhohat,numel(p));

% Leave-one-dataset-out cross-validation tests generalization to unseen datasets.
if doLODO
    cv = leave_one_dataset_out_cv_fast_fixed_hcrit( ...
        rho_fit,T_fit,h_fit,dataset_fit,p, ...
        h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid, ...
        Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm);
else
    cv.rmse = NaN;
    cv.r2 = NaN;
end

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
set(fig1,'Units','inches','Position',[1.1 1.0 10.5 8.8]);
tl = tiledlayout(fig1,2,2,'TileSpacing','compact','Padding','compact');

ax0_tmp = nexttile(tl,1);
drawnow

set(ax0_tmp,'Units','inches');
tilePos = ax0_tmp.Position;
delete(ax0_tmp)

img_map = imread(fullfile(figureDir, ...
    'Map_density_lab.png'));
[imgH,imgW,~] = size(img_map);
dpi = 300;
imgW_in = imgW / dpi;
imgH_in = imgH / dpi;

dx = -0.18;
dy = +0.01;
left = tilePos(1) + dx;
top = tilePos(2) + tilePos(4) + dy;

ax0 = axes(fig1, ...
    'Units','inches', ...
    'Position',[left, top-imgH_in, imgW_in, imgH_in]);

image(ax0,img_map)
axis(ax0,'image')
axis(ax0,'off')
set(ax0,'YDir','reverse')

ax1 = nexttile(tl,2);
hold(ax1,'on')
scatter(ax1,t(is_fyi_all),rho(is_fyi_all),38,T(is_fyi_all),'o','filled','MarkerEdgeColor','flat')
scatter(ax1,t(is_old_all),rho(is_old_all),38,T(is_old_all),'o','MarkerFaceColor','none','MarkerEdgeColor','flat','LineWidth',1.0)
scatter(ax1,t(is_sva_all),rho(is_sva_all),38,T(is_sva_all),'>','filled','MarkerEdgeColor','flat')

idxGuide = ice_type == "FYI" & ~is_svalbard & ~isnan(rho) & ~isnat(t);
tGuide = t(idxGuide);
rhoGuide = rho(idxGuide);
t_start_plot = datetime(2020,1,1);
t_end_plot = datetime(2021,1,1);
tRef = datetime(2020,1,1);
xGuide = days(tGuide - tRef);
idxDec = month(tGuide) == 12;
idxJan = month(tGuide) == 1;
xWrap = [xGuide(idxDec)-366; xGuide; xGuide(idxJan)+366];
rhoWrap = [rhoGuide(idxDec); rhoGuide; rhoGuide(idxJan)];

% Robust LOESS smoothing used only for visualization of the FYI seasonal density cycle; it is not used in model fitting.
if numel(xWrap) >= 5
    [xU,~,ic] = unique(xWrap);
    rhoDay = accumarray(ic,rhoWrap,[],@median);
    [xU,isrt] = sort(xU);
    rhoDay = rhoDay(isrt);
    tCurve = (t_start_plot:days(1):t_end_plot-days(1)).';
    xCurve = days(tCurve - tRef);
    rhoSmooth = smooth(xU,rhoDay,0.25,'rloess');
    rhoCurve = interp1(xU,rhoSmooth,xCurve,'pchip',nan);
    p_mean = plot(ax1,tCurve,rhoCurve,'--','Color',[0.15 0.15 0.15],'LineWidth',1.2);
else
    p_mean = plot(ax1,nan,nan,'--','Color',[0.15 0.15 0.15],'LineWidth',1.2);
end

p_fyi_ca = plot(ax1,nan,nan,'o','MarkerSize',ms_ca_fyi,'MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','none');
p_syi_ca = plot(ax1,nan,nan,'o','MarkerSize',ms_ca_syi,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',lw_open,'LineStyle','none');
p_sva = plot(ax1,nan,nan,'>','MarkerSize',ms_sva,'MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','none');
legend(ax1,[p_fyi_ca p_syi_ca p_sva p_mean],{'FYI','SYI & MYI','Svalbard','FYI seasonal mean'},'Location','southwest','Box','off','FontSize',fs_txt,'Interpreter','tex')
colormap(ax1,lipari)
clim(ax1,[Tmin_time Tmax_time])
cb1 = colorbar(ax1);
ylabel(cb1, 'Sea-ice temperature (°C)');
xlabel(ax1,'Month')
ylabel(ax1,'Sea-ice density \rho (kg m^{-3})','Interpreter','tex')
set(ax1,'FontSize',fs_ax,'FontWeight','normal')
xlim(ax1,[t_start_plot t_end_plot])
datetick(ax1, 'x', 'mmm', 'keepticks');
ylim(ax1,[850 920])
box(ax1,'on')

fsz = 10;
label_dx = -105;
label_dy = -3.87;
text(ax0, ...
    label_dx, imgH + label_dy*dpi, '(a)', ...
    'Units','data', ...
    'FontSize',fsz, ...
    'FontWeight','normal', ...
    'Clipping','off');

text(ax1,-0.145,1.02,'(b)','Units','normalized','FontSize',10,'FontWeight','normal')

ax2 = nexttile(tl,3);
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

leg2.Units = 'normalized';
leg2.ItemTokenSize = [14 4];
leg2.Position = [0.015 0.23 0.23 0.13];

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
xlabel(ax2,'Sea-ice temperature{\it T_i} (°C)')
ylabel(ax2,'Sea-ice density \rho (kg m^{-3})','Interpreter','tex')
set(ax2,'FontSize',fs_ax,'FontWeight','normal')
xlim(ax2,[-10 0])
ylim(ax2,[850 920])
box(ax2,'on')
text(ax2,-0.145,1.02,'(c)','Units','normalized','FontSize',10,'FontWeight','normal')

ax3 = nexttile(tl,4);
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

p_fyi3 = plot(ax3,nan,nan,'o', ...
    'MarkerSize',ms_ca_fyi, ...
    'MarkerFaceColor','k', ...
    'MarkerEdgeColor','k', ...
    'LineStyle','none');

p_syi3 = plot(ax3,nan,nan,'o', ...
    'MarkerSize',ms_ca_syi, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor','k', ...
    'LineWidth',lw_open, ...
    'LineStyle','none');

p_sva3 = plot(ax3,nan,nan,'>', ...
    'MarkerSize',ms_sva, ...
    'MarkerFaceColor','k', ...
    'MarkerEdgeColor','k', ...
    'LineStyle','none');

leg3 = legend(ax3, ...
    [p_fyi3 p_syi3 p_sva3 p_T1 p_T2 p_T3 h_ci3], ...
    {'FYI','SYI \& MYI','Svalbard', ...
     sprintf('$T_i = %.1f^{\\circ}\\mathrm{C}$',T_sel(1)), ...
     sprintf('$T_i = %.1f^{\\circ}\\mathrm{C}$',T_sel(2)), ...
     sprintf('$T_i = %.1f^{\\circ}\\mathrm{C}$',T_sel(3)), ...
     '95\% bootstrap CI'}, ...
    'Location','southeast', ...
    'Box','off', ...
    'FontSize',fs_txt);

leg3.AutoUpdate = 'off';
set(leg3,'FontSize',fs_txt,'Interpreter','latex')
xlabel(ax3,'Sea-ice thickness{\it h_i} (m)')
ylabel(ax3,'Sea-ice density \rho (kg m^{-3})','Interpreter','tex')
set(ax3,'FontSize',fs_ax,'FontWeight','normal')
xlim(ax3,[0 3.5])
ylim(ax3,[850 920])
box(ax3,'on')
text(ax3,-0.145,1.02,'(d)','Units','normalized','FontSize',10,'FontWeight','normal')

outputFigure = fullfile(figureDir,'Fig1.png');

set(findall(fig1,'Type','axes'),'Toolbar',[])

exportgraphics(fig1, outputFigure,'Resolution',300)

close(fig1)

fprintf('Generated Figure 1 using %s density (%d observations).\n', rho_names{densityMode}, N)
fprintf('Saved figure to:\n%s\n', outputFigure)

%% Helpers
function fit = fit_fast_fixed_hcrit_model(rho,T,h,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid,Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm)
% Grid-search optimization of the temperature-thickness density model.

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
% Evaluate the fitted density parameterization.

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

% Piecewise-linear cold-regime density parameterization.
% Density increases to a maximum at h_cold_peak and then decreases
% gradually toward h_cold_high.
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
% Leave-one-dataset-out cross-validation treating each campaign as an independent observational unit.

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
% Compute standard regression performance metrics.

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
% Dataset-level bootstrap estimation of parameter uncertainty.

n = numel(rho);
B = nan(nboot,5);
groups = unique(dataset_id);
G = numel(groups);

for b = 1:nboot
    % Resample entire observational datasets rather than individual observations to preserve within-campaign dependence structures.
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
% Plot 95% bootstrap confidence intervals for model predictions.

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

end
