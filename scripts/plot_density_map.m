function plot_density_map(repoRoot)

close all

finalDir = fullfile(repoRoot, 'data', 'final');
figureDir = fullfile(repoRoot, 'figures');

if ~exist(figureDir, 'dir')
    mkdir(figureDir)
end

inputFile = fullfile(finalDir, ...
    'snow_model_matchups_with_models.mat');

load(inputFile,'D')

if ~exist('D','var') || ~istable(D)
    error('D not found in snow_model_matchups.mat')
end

Dmap = table();

Dmap.rho = double(D.density_lab_kgm3);
Dmap.lat = double(D.lat);
Dmap.lon = double(D.lon);

idx = ~isnan(Dmap.lat) & ~isnan(Dmap.lon) & ~isnan(Dmap.rho);
Dmap = Dmap(idx,:);

if isempty(Dmap)
    error('No valid rows with lat, lon, and density.')
end

rho_min = 870;
rho_max = 920;

ncol_hi = 256;
colormapDir = fullfile(repoRoot, ...
    'data', 'colormaps');
S = load(fullfile(colormapDir,'buda.mat'));
c_hi = flipud(S.buda);

rho_scaled = (Dmap.rho - rho_min) ./ (rho_max - rho_min);
hi_idx = round(1 + rho_scaled * (ncol_hi - 1));
hi_idx = max(1, min(ncol_hi, hi_idx));
hi_rgb = c_hi(hi_idx, :);

lon0 = 20;
lat0 = min(89, max(Dmap.lat, [], 'omitnan') - 0.5);

lon_span = max(Dmap.lon) - min(Dmap.lon);
lat_span = max(Dmap.lat) - min(Dmap.lat);

radius = max(6.2, 0.36 * max(lon_span/2, lat_span*2.0));
radius = min(radius, 16);
radius = radius * 1.02;
fs = 9.5;

fig = figure;
tiledlayout(1,1,'TileSpacing','compact','Padding','compact')
nexttile

m_proj('stereographic', ...
    'lat', lat0, ...
    'long', lon0, ...
    'radius', radius, ...
    'rectbox', 'on');

depth_levels = -6000:200:0;
[~,~] = m_etopo2('contourf', depth_levels, 'edgecolor', 'none');
shading flat
% Sdeep = load(fullfile(colormapDir,'deep.mat'));
% colormap(flipud(Sdeep.deep))
Sgray = load(fullfile(colormapDir,'grayC.mat'));
colormap(Sgray.grayC)
clim([-6000 0])

hold on

% Land patch without edge, then coastline on top
m_gshhs_i('patch', [0.78 0.78 0.78], 'edgecolor', 'none');
m_gshhs_i('color', [0.5 0.5 0.5], 'linewi', 0.5);

% Sampling markers
h_sc = m_scatter(Dmap.lon, Dmap.lat, 58-35, hi_rgb, ...
    'filled', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.2);

m_grid('xtick', [-150 -100 -50 0 50 100], ...
       'ytick', [78 80 82 84 86 88], ...
       'tickdir', 'in', ...
       'linestyle', ':', ...
       'color', [0 0 0], ...
       'linewidth', 0.01, ...
       'fontsize', fs, ...
       'box', 'on', ...
       'linewi', 0.3);

drawnow

hx = findall(gca,'tag','m_grid_xticklabel');
hy = findall(gca,'tag','m_grid_yticklabel');

if ~isempty(hx)
    px = get(hx,'position');
    if ~iscell(px), px = {px}; end
    yx = cellfun(@(p) p(2), px);
    ytop_x = max(yx);
    set(hx(yx > ytop_x - 0.02*range(yx)), 'String', '');
end

if ~isempty(hy)
    py = get(hy,'position');
    if ~iscell(py), py = {py}; end
    yy = cellfun(@(p) p(2), py);
    ytop_y = max(yy);
    set(hy(yy > ytop_y - 0.02*range(yy)), 'String', '');
end

% --- Basin labels ---
m_text(60, 85.0, 'Nansen Basin', ...
    'fontsize', fs-1.5, ...
    'color', [0.65 0.65 0.65], ...
    'fontangle', 'italic', ...
    'horizontalalignment', 'center', ...
    'rotation', 35);
m_text(10, 87.9+0.3, 'Amundsen Basin', ...
    'fontsize', fs-1.5, ...
    'color', [0.65 0.65 0.65], ...
    'fontangle', 'italic', ...
    'horizontalalignment', 'center', ...
    'rotation', 35);
m_text(43+4, 79, 'Barents Sea', ...
    'fontsize', fs-1.5, ...
    'color', [0.65 0.65 0.65], ...
    'fontangle', 'italic', ...
    'horizontalalignment', 'center', ...
    'rotation', 35);
m_text(-8.5, 80.0, 'Fram Strait', ...
    'fontsize', fs-1.5, ...
    'color', [0.4 0.4 0.4], ...
    'fontangle', 'italic', ...
    'horizontalalignment', 'center', ...
    'rotation', 70);

% Bring markers back to top
uistack(h_sc, 'top');

% Make tick-label text black
set(findall(gca,'tag','m_grid_xticklabel'),'Color','k');
set(findall(gca,'tag','m_grid_yticklabel'),'Color','k');

% -----------------------------
% Depth colorbar on right
% -----------------------------
[ax_depth,~] = m_contfbar(1.05, [0.25 0.82], [-6000 0], depth_levels, ...
    'edgecolor', 'none', ...
    'endpiece', 'no', ...
    'axfrac', 0.022, ...
    'tickdir', 'out');

t_depth = title(ax_depth, 'Depth (km)', 'fontsize', fs-1, 'fontweight', 'normal');
t_depth.Units = 'normalized'; t_depth.Position(2) = 1.05;
set(ax_depth, ...
    'YTick', [-6000 -4000 -2000 0], ...
    'YTickLabel', {'6','4','2','0'}, ...
    'FontSize', fs-1);

% -----------------------------
% Ice-density colorbar on left
% -----------------------------
img_hi = permute(c_hi, [1 3 2]);
ax_hi = axes;
image(ax_hi, img_hi);

set(ax_hi, ...
    'YDir', 'normal', ...
    'YAxisLocation', 'right', ...
    'XTick', [], ...
    'TickDir', 'out', ...
    'FontSize', fs, ...
    'Box', 'on', ...
    'LineWidth', 0.8);

% rho_ticks = linspace(rho_min, rho_max, 4);
rho_ticks = 870:10:920;
tick_idx = linspace(1, ncol_hi, numel(rho_ticks));

set(ax_hi, ...
    'YTick', tick_idx, ...
    'YTickLabel', compose('%.0f', rho_ticks));

t_rho = title(ax_hi, 'Sea-ice density (kg m^{-3})', 'FontSize', fs, 'FontWeight', 'normal');
t_rho.Units = 'normalized'; t_rho.Position(2) = 1.05;

% MANUAL COLORBAR POSITIONS: [left bottom width height]
cb_bottom = 0.68;
cb_height = 0.22;
cb_width  = 0.018;

left_cb_left  = 0.75-0.35;
right_cb_left = 0.75+0.08;

set(ax_hi,    'Position', [left_cb_left,  cb_bottom, cb_width, cb_height]);
set(ax_depth, 'Position', [right_cb_left, cb_bottom, cb_width, cb_height]);

set(fig, 'Units', 'inches', 'Position', [3 5 4.39 4.06]);

outputFigure = fullfile(figureDir, ...
    'Map_density_lab.png');

set(findall(fig,'Type','axes'),'Toolbar',[])

exportgraphics(fig, outputFigure, 'Resolution', 300);

close(fig)

fprintf('Saved density map to:\n%s\n', ...
    outputFigure)

end
