clear
close all
clc

repoRoot = 'C:\Users\evsalg001\Documents\MATLAB\Density parametrization\script';

finalDir = fullfile(repoRoot, ...
    'data', 'final');

colormapDir = fullfile(repoRoot, ...
    'data', 'colormaps');

inputFile = fullfile(finalDir, ...
    'snow_model_matchups_with_models.mat');

load(inputFile,'D')

S = load(fullfile(colormapDir,'lipari.mat'));
lipari = S.lipari;

densityMode = 1; % 1 = lab, 2 = connected, 3 = disconnected
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

if ~isdatetime(D.date)
    D.date = datetime(D.date);
end

% Freeze-up definition for this test:
% FYI cores sampled in October or November.
ice_type_all = upper(string(D.ice_age));
is_freezeup_all = ice_type_all == "FYI" & ismember(month(D.date),[10 11]);

% Keep these only as additional flags. They are not used for the freeze-up
% definition above.
dataset_all = string(D.dataset);
is_svalbard_all = contains(lower(dataset_all),"svalbard");

t = D.date;
t.Year = 2020;
T = D.temperature_C;
rho = D.(rho_col);
h = D.ice_thickness_m;
T = D.temperature_C;
dataset_id = string(D.dataset);
ice_type = upper(string(D.ice_age));
is_freezeup = is_freezeup_all;
is_svalbard = is_svalbard_all;

valid = ~isnat(t) & isfinite(T) & isfinite(rho) & isfinite(h) & T < 0 & h > 0;

t = t(valid);
T = T(valid);
rho = rho(valid);
h = h(valid);
dataset_id = dataset_id(valid);
ice_type = ice_type(valid);
is_freezeup = is_freezeup(valid);
is_svalbard = is_svalbard(valid);

% Plot groups. Freeze-up cores are removed from the regular FYI group so
% they can be drawn on top with a special marker.
is_fyi = ice_type == "FYI" & ~is_svalbard & ~is_freezeup;
is_old = ice_type ~= "FYI" & ~is_svalbard & ~is_freezeup;
is_sva = is_svalbard & ~is_freezeup;

% Fit exactly as in the script-2 approach, but now including freeze-up cores
% in the all-data fit. Svalbard is excluded from the fit, matching the logic
% used for the main non-Svalbard model comparison.
fit_mask_all = ~is_svalbard;

T_fit = T(fit_mask_all);
rho_fit = rho(fit_mask_all);
h_fit = h(fit_mask_all);
dataset_fit = dataset_id(fit_mask_all);
is_freezeup_fit = is_freezeup(fit_mask_all);
N_all = numel(rho_fit);

% Cold branch for T_i <= Tcrit.
h_cold_peak = 1.10;
rho_cold_peak = 911;
h_cold_high = 3.50;
rho_cold_high = 908.5;
rho_cold_0_grid = 905:1:912;

% Warm branch for T_i > Tcrit.
Tcrit_grid = -2.8:0.1:-2.1;
hcrit_warm_grid = 1.8:0.1:2.3;
rho_warm_plateau_grid = 895:1:905;

min_n_cold = 8;
min_n_warm = 8;
nboot = 100;
rng(1)

fit = fit_fast_fixed_hcrit_model( ...
    rho_fit,T_fit,h_fit, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid, ...
    Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm);

p = [fit.Tcrit fit.rho_warm_plateau fit.warm_slope fit.hcrit_warm fit.rho_cold_0];

rhohat_all = predict_fast_fixed_hcrit_model( ...
    T_fit,h_fit,p,h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high);

metrics_all = regression_metrics(rho_fit,rhohat_all,numel(p));

idx_eval_no_freezeup = ~is_freezeup_fit;
metrics_no_freezeup = regression_metrics( ...
    rho_fit(idx_eval_no_freezeup),rhohat_all(idx_eval_no_freezeup),numel(p));

N_no_freezeup = sum(idx_eval_no_freezeup);
N_freezeup = sum(is_freezeup_fit);

h_freezeup_fit = h_fit(is_freezeup_fit);
T_freezeup_fit = T_fit(is_freezeup_fit);

h_freezeup_p95 = prctile(h_freezeup_fit,95);
T_freezeup_p05 = prctile(T_freezeup_fit,5);

fprintf('\nFreeze-up core physical range in fit domain:\n')
fprintf('Thickness median = %.2f m, range = %.2f--%.2f m, 95th percentile = %.2f m\n', ...
    median(h_freezeup_fit,'omitnan'), ...
    min(h_freezeup_fit,[],'omitnan'), ...
    max(h_freezeup_fit,[],'omitnan'), ...
    h_freezeup_p95)

fprintf('Temperature median = %.1f C, range = %.1f--%.1f C, 5th percentile = %.1f C\n', ...
    median(T_freezeup_fit,'omitnan'), ...
    min(T_freezeup_fit,[],'omitnan'), ...
    max(T_freezeup_fit,[],'omitnan'), ...
    T_freezeup_p05)

B_fit = bootstrap_fast_fixed_hcrit_model( ...
    rho_fit,T_fit,h_fit,dataset_fit,nboot, ...
    h_cold_peak,rho_cold_peak,h_cold_high,rho_cold_high,rho_cold_0_grid, ...
    Tcrit_grid,rho_warm_plateau_grid,hcrit_warm_grid,min_n_cold,min_n_warm);

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
ms_freezeup = 7.2;
lw_open = 1.0;
Tmin_time = -8;
Tmax_time = 0;

fig1 = figure;
set(fig1,'Units','inches','Position',[1.1 1.0 10.5 4.4]);
tl = tiledlayout(fig1,1,2,'TileSpacing','compact','Padding','compact');

% Panel (b): density vs month, with freeze-up cores highlighted.
ax1 = nexttile(tl,1);
hold(ax1,'on')
scatter(ax1,t(is_fyi),rho(is_fyi),38,T(is_fyi),'o','filled','MarkerEdgeColor','flat')
scatter(ax1,t(is_old),rho(is_old),38,T(is_old),'o','MarkerFaceColor','none','MarkerEdgeColor','flat','LineWidth',1.0)
scatter(ax1,t(is_sva),rho(is_sva),38,T(is_sva),'>','filled','MarkerEdgeColor','flat')
scatter(ax1,t(is_freezeup),rho(is_freezeup),64,T(is_freezeup),'d','filled','MarkerEdgeColor','k','LineWidth',0.9)

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
p_freezeup = plot(ax1,nan,nan,'d','MarkerSize',ms_freezeup,'MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','none');
legend(ax1,[p_fyi_ca p_syi_ca p_sva p_freezeup p_mean], ...
    {'FYI','SYI & MYI','Svalbard','Freeze-up FYI, Oct-Nov','FYI seasonal mean'}, ...
    'Location','southwest','Box','off','FontSize',fs_txt,'Interpreter','tex')
colormap(ax1,lipari)
clim(ax1,[Tmin_time Tmax_time])
cb1 = colorbar(ax1);
ylabel(cb1, 'Sea-ice temperature (°C)');
xlabel(ax1,'Month')
ylabel(ax1,'Sea-ice density \rho (kg m^{-3})','Interpreter','tex')
set(ax1,'FontSize',fs_ax,'FontWeight','normal')
xlim(ax1,[t_start_plot t_end_plot])
datetick(ax1,'x','mmm','keepticks')
ylim(ax1,[850 920])
box(ax1,'on')
text(ax1,-0.145,1.02,'(b)','Units','normalized','FontSize',10,'FontWeight','normal')

% Panel (d): density vs thickness, with same all-data fit but skill also
% evaluated after excluding freeze-up cores.
ax3 = nexttile(tl,2);
hold(ax3,'on')
scatter(ax3,h(is_fyi),rho(is_fyi),ms,T(is_fyi),'filled','MarkerEdgeColor','flat','LineWidth',0.35)
scatter(ax3,h(is_old),rho(is_old),ms,T(is_old),'MarkerFaceColor','none','MarkerEdgeColor','flat','LineWidth',1.0)
scatter(ax3,h(is_sva),rho(is_sva),ms,T(is_sva),'>','filled','MarkerEdgeColor','flat','LineWidth',0.35)
scatter(ax3,h(is_freezeup),rho(is_freezeup),64,T(is_freezeup),'d','filled','MarkerEdgeColor','k','LineWidth',0.9)
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

p_freezeup3 = plot(ax3,nan,nan,'d', ...
    'MarkerSize',ms_freezeup, ...
    'MarkerFaceColor','k', ...
    'MarkerEdgeColor','k', ...
    'LineStyle','none');

leg3 = legend(ax3, ...
    [p_fyi3 p_syi3 p_sva3 p_freezeup3 p_T1 p_T2 p_T3 h_ci3], ...
    {'FYI','SYI \& MYI','Svalbard','Freeze-up FYI, Oct-Nov', ...
     sprintf('$T_i = %.1f^{\\circ}\\mathrm{C}$',T_sel(1)), ...
     sprintf('$T_i = %.1f^{\\circ}\\mathrm{C}$',T_sel(2)), ...
     sprintf('$T_i = %.1f^{\\circ}\\mathrm{C}$',T_sel(3)), ...
     '95\% bootstrap CI'}, ...
    'Location','southeast', ...
    'Box','off', ...
    'FontSize',fs_txt);

leg3.AutoUpdate = 'off';
set(leg3,'FontSize',fs_txt,'Interpreter','latex')

skillTxt = sprintf([ ...
    'Same fit evaluated two ways' newline ...
    'All data: R^2 = %.2f, RMSE = %.1f kg m^{-3}, N = %d' newline ...
    'No freeze-up: R^2 = %.2f, RMSE = %.1f kg m^{-3}, N = %d' newline ...
    'Improvement: dR^2 = %.2f, dRMSE = %.1f kg m^{-3}'], ...
    metrics_all.r2,metrics_all.rmse,N_all, ...
    metrics_no_freezeup.r2,metrics_no_freezeup.rmse,N_no_freezeup, ...
    metrics_no_freezeup.r2 - metrics_all.r2, ...
    metrics_all.rmse - metrics_no_freezeup.rmse);

text(ax3,0.02,0.01,skillTxt, ...
    'Units','normalized', ...
    'Interpreter','tex', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'FontSize',fs_txt + 0.8)

xlabel(ax3,'Sea-ice thickness{\it h_i} (m)')
ylabel(ax3,'Sea-ice density \rho (kg m^{-3})','Interpreter','tex')
set(ax3,'FontSize',fs_ax,'FontWeight','normal')
xlim(ax3,[0 3.5])
ylim(ax3,[850 920])
box(ax3,'on')
text(ax3,-0.145,1.02,'(d)','Units','normalized','FontSize',10,'FontWeight','normal')

set(findall(fig1,'Type','axes'),'Toolbar',[])

fprintf('Generated freeze-up skill test using %s density.\n',rho_names{densityMode})
fprintf('Freeze-up definition: FYI cores in October or November.\n')
fprintf('All-data fit:       R2 = %.3f, RMSE = %.2f kg m-3, N = %d\n', ...
    metrics_all.r2,metrics_all.rmse,N_all)
fprintf('No freeze-up eval:  R2 = %.3f, RMSE = %.2f kg m-3, N = %d\n', ...
    metrics_no_freezeup.r2,metrics_no_freezeup.rmse,N_no_freezeup)
fprintf('Freeze-up cores in fit domain: N = %d\n',N_freezeup)
fprintf('Skill improvement after removing freeze-up from evaluation: dR2 = %.3f, dRMSE = %.2f kg m-3\n', ...
    metrics_no_freezeup.r2 - metrics_all.r2,metrics_all.rmse - metrics_no_freezeup.rmse)
fprintf('Freeze-up physical support: h_95 = %.2f m, T_05 = %.1f C\n', ...
    h_freezeup_p95,T_freezeup_p05)

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

rhoCold = nan(size(h));

idx1 = h <= h_cold_peak;
rhoCold(idx1) = rho_cold_0 + ...
    (rho_cold_peak - rho_cold_0).*h(idx1)./h_cold_peak;

idx2 = h > h_cold_peak;
rhoCold(idx2) = rho_cold_peak + ...
    (rho_cold_high - rho_cold_peak).* ...
    (min(h(idx2),h_cold_high) - h_cold_peak)./(h_cold_high - h_cold_peak);
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
