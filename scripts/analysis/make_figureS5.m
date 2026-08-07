function make_figureS5(repoRoot)
%
% MAKE_FIGURES5
%
% Generates Supplementary Figure S5: CryoSat‑2 sea‑ice thickness differences
% between This Study and two reference density parameterizations:
%
%   • A2010 – Alexandrov et al. (2010)
%   • J2022 – Jutila et al. (2022)
%
% Figure layout:
%   Row 1: February 2011–2023 climatology
%   Row 2: July 2011–2023 climatology
%
%   Columns 1–2: This Study minus A2010 (map + PDF)
%   Columns 3–4: This Study minus J2022 (map + PDF)
%
% Interpretation:
%   Positive values → This Study produces thicker ice
%   Negative values → This Study produces thinner ice
%
% Required input files (placed in data/cs2/):
%   uit_sit_v2northpolarstereo_80km_feb_2011-2023_v3p0.nc
%   uit_sit_v2northpolarstereo_80km_jul_2011-2023_v3p0.nc
%
% Output:
%   figures/FigS5.png
%

% Suppress Mapping Toolbox warnings
warning('off','map:removing:mfwdtran')
warning('off','map:removing:minvtran')

% Load February CryoSat‑2 SIT and SIC fields.
filename_feb = fullfile(repoRoot, 'data', 'cs2', ...
    'uit_sit_v2northpolarstereo_80km_feb_2011-2023_v3p0.nc');

grid_lat = ncread(filename_feb,'latitude');
grid_lon = ncread(filename_feb,'longitude');
sic_feb = ncread(filename_feb,'sic');

sit_A2010_feb = ncread(filename_feb,'sit_A2010');
sit_J2021_feb = ncread(filename_feb,'sit_J2021');
sit_S2026_feb = ncread(filename_feb,'sit_S2026');

% Load July CryoSat‑2 SIT and SIC fields
filename_jul = fullfile(repoRoot, 'data', 'cs2', ...
    'uit_sit_v2northpolarstereo_80km_jul_2011-2023_v3p0.nc');

sic_jul = ncread(filename_jul,'sic');

sit_A2010_jul = ncread(filename_jul,'sit_A2010');
sit_J2021_jul = ncread(filename_jul,'sit_J2021');
sit_S2026_jul = ncread(filename_jul,'sit_S2026');

% ColorBrewer RdBu diverging colormap.
rdBu11 = [ ...
    103   0  31
    178  24  43
    214  96  77
    244 165 130
    253 219 199
    247 247 247
    209 229 240
    146 197 222
     67 147 195
     33 102 172
      5  48  97] / 255;

% Interpolate the 11-color reference table to a smooth 64-color map.
x = linspace(0,1,size(rdBu11,1));
xq = linspace(0,1,64);
cmap = interp1(x,rdBu11,xq,'linear');

% Define thickness difference range for maps and histograms

sit_diff_int = -1:0.02:1;

% Create 2×4 tiled figure layout
% Row 1: February differences
% Row 2: July differences
% Columns:
%   1–2 → This Study − A2010 (map + PDF)
%   3–4 → This Study − J2022 (map + PDF)
figS5 = figure(1);
clf(figS5)

set(figS5, ...
    'Units','pixels', ...
    'Position',[100 100 1617 768], ...
    'Color','white')

tl = tiledlayout(figS5,2,4, ...
    'TileSpacing','compact');

% (a) February: This Study minus A2010 (map)
ax1 = nexttile(tl);

% evalc captures command-window output generated internally by ncpolarm,
% including the repeated "Northern Hemisphere" message.
evalc("ncpolarm('lat',68,'lon',0);");

temp = nanmean(sit_S2026_feb - sit_A2010_feb,3);

% Fill missing values near the pole where CryoSat‑2 has no coverage.
temp2 = inpaint_nans(temp,1);
idxPole = grid_lat > 88 & isnan(temp);
temp(idxPole) = temp2(idxPole);

% Apply spatial and sea-ice-concentration masks.
temp(grid_lat > 88.5 | ...
     temp == 0 | ...
     nanmean(sic_feb,3) < 15) = NaN;

contourfm(grid_lat,grid_lon,temp,sit_diff_int, ...
    'LineStyle','none')

contourm(grid_lat,grid_lon,nanmean(sic_feb,3), ...
    [15 15], ...
    'k', ...
    'Fill','off', ...
    'LineWidth',0.5)

clim(ax1,[min(sit_diff_int)-0.01 max(sit_diff_int)])
colormap(ax1,cmap)

title(ax1,'This Study minus A2010', ...
    'FontSize',16, ...
    'FontWeight','normal')

text(ax1,-0.05,0.5,'Feb 2011-2023', ...
    'Units','normalized', ...
    'Rotation',90, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontSize',16, ...
    'FontWeight','normal')

text(ax1,-0.13,1.02,'(a)', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','normal', ...
    'Clipping','off')

% (b) February: Thickness difference distribution (PDF)
ax2 = nexttile(tl);
cla(ax2)
box(ax2,'on')
grid(ax2,'on')
hold(ax2,'on')

xlabel(ax2,'Sea-ice thickness difference (m)')
ylabel(ax2,'PDF')

histogram(ax2, ...
    sit_S2026_feb - sit_A2010_feb, ...
    sit_diff_int, ...
    'Normalization','pdf', ...
    'LineStyle','none', ...
    'FaceColor',cmap(end,:), ...
    'FaceAlpha',0.9)

xlim(ax2,[-1 1])
set(ax2,'FontSize',12)

text(ax2,-0.21,1.02,'(b)', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','normal', ...
    'Clipping','off')

% February: This Study minus J2022 map.
ax3 = nexttile(tl);
evalc("ncpolarm('lat',68,'lon',0);");

temp = nanmean(sit_S2026_feb - sit_J2021_feb,3);

temp2 = inpaint_nans(temp,1);
idxPole = grid_lat > 88 & isnan(temp);
temp(idxPole) = temp2(idxPole);

% Mask unreliable grid cells:
%   • lat > 88.5° → outside reliable CS2 coverage
%   • temp == 0 → missing SIT values
%   • SIC < 15% → standard threshold for sea‑ice presence

temp(grid_lat > 88.5 | ...
     temp == 0 | ...
     nanmean(sic_feb,3) < 15) = NaN;

contourfm(grid_lat,grid_lon,temp,sit_diff_int, ...
    'LineStyle','none')

contourm(grid_lat,grid_lon,nanmean(sic_feb,3), ...
    [15 15], ...
    'k', ...
    'Fill','off', ...
    'LineWidth',0.5)

clim(ax3,[min(sit_diff_int)-0.01 max(sit_diff_int)])
colormap(ax3,cmap)

title(ax3,'This Study minus J2022', ...
    'FontSize',16, ...
    'FontWeight','normal')

text(ax3,-0.13,1.02,'(c)', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','normal', ...
    'Clipping','off')

% February: distribution of This Study minus J2022.
ax4 = nexttile(tl);
cla(ax4)
box(ax4,'on')
grid(ax4,'on')
hold(ax4,'on')

xlabel(ax4,'Sea-ice thickness difference (m)')
ylabel(ax4,'PDF')

histogram(ax4, ...
    sit_S2026_feb - sit_J2021_feb, ...
    sit_diff_int, ...
    'Normalization','pdf', ...
    'LineStyle','none', ...
    'FaceColor',cmap(end,:), ...
    'FaceAlpha',0.9)

xlim(ax4,[-1 1])
set(ax4,'FontSize',12)

text(ax4,-0.23,1.02,'(d)', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','normal', ...
    'Clipping','off')

% (e) July: This Study minus A2010 (map)
ax5 = nexttile(tl);
evalc("ncpolarm('lat',68,'lon',0);");

temp = nanmean(sit_S2026_jul - sit_A2010_jul,3);

temp2 = inpaint_nans(temp,1);
idxPole = grid_lat > 88 & isnan(temp);
temp(idxPole) = temp2(idxPole);

temp(grid_lat > 88.5 | ...
     temp == 0 | ...
     nanmean(sic_jul,3) < 15) = NaN;

contourfm(grid_lat,grid_lon,temp,sit_diff_int, ...
    'LineStyle','none')

contourm(grid_lat,grid_lon,nanmean(sic_jul,3), ...
    [15 15], ...
    'k', ...
    'Fill','off', ...
    'LineWidth',0.5)

clim(ax5,[min(sit_diff_int)-0.01 max(sit_diff_int)])
colormap(ax5,cmap)

title(ax5,'This Study minus A2010', ...
    'FontSize',16, ...
    'FontWeight','normal')

text(ax5,-0.05,0.5,'July 2011-2023', ...
    'Units','normalized', ...
    'Rotation',90, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontSize',16, ...
    'FontWeight','normal')

text(ax5,-0.13,1.02,'(e)', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','normal', ...
    'Clipping','off')

% July: distribution of This Study minus A2010.
ax6 = nexttile(tl);
cla(ax6)
box(ax6,'on')
grid(ax6,'on')
hold(ax6,'on')

xlabel(ax6,'Sea-ice thickness difference (m)')
ylabel(ax6,'PDF')

histogram(ax6, ...
    sit_S2026_jul - sit_A2010_jul, ...
    sit_diff_int, ...
    'Normalization','pdf', ...
    'LineStyle','none', ...
    'FaceColor',cmap(end,:), ...
    'FaceAlpha',0.9)

xlim(ax6,[-1 1])
set(ax6,'FontSize',12)

text(ax6,-0.21,1.02,'(f)', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','normal', ...
    'Clipping','off')

% July: This Study minus J2022 map.
ax7 = nexttile(tl);
evalc("ncpolarm('lat',68,'lon',0);");

temp = nanmean(sit_S2026_jul - sit_J2021_jul,3);

temp2 = inpaint_nans(temp,1);
idxPole = grid_lat > 88 & isnan(temp);
temp(idxPole) = temp2(idxPole);

temp(grid_lat > 88.5 | ...
     temp == 0 | ...
     nanmean(sic_jul,3) < 15) = NaN;

contourfm(grid_lat,grid_lon,temp,sit_diff_int, ...
    'LineStyle','none')

contourm(grid_lat,grid_lon,nanmean(sic_jul,3), ...
    [15 15], ...
    'k', ...
    'Fill','off', ...
    'LineWidth',0.5)

clim(ax7,[min(sit_diff_int)-0.01 max(sit_diff_int)])
colormap(ax7,cmap)

title(ax7,'This Study minus J2022', ...
    'FontSize',16, ...
    'FontWeight','normal')

c = colorbar(ax7,'south', ...
    'FontSize',13);

c.AxisLocation = 'out';
c.Position(2) = c.Position(2) - 0.04;
c.Label.String = 'Sea-ice thickness difference (m)';

text(ax7,-0.13,1.02,'(g)', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','normal', ...
    'Clipping','off')

% July: distribution of This Study minus J2022.
ax8 = nexttile(tl);
cla(ax8)
box(ax8,'on')
grid(ax8,'on')
hold(ax8,'on')

xlabel(ax8,'Sea-ice thickness difference (m)')
ylabel(ax8,'PDF')

histogram(ax8, ...
    sit_S2026_jul - sit_J2021_jul, ...
    sit_diff_int, ...
    'Normalization','pdf', ...
    'LineStyle','none', ...
    'FaceColor',cmap(end,:), ...
    'FaceAlpha',0.9)

xlim(ax8,[-1 1])
set(ax8,'FontSize',12)

text(ax8,-0.23,1.02,'(h)', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','normal', ...
    'Clipping','off')

drawnow

% Export Figure S5 as a 300 dpi PNG.
outputFile = fullfile(repoRoot, 'figures', 'FigS5.png');
exportgraphics(figS5, outputFile, 'Resolution',300, 'BackgroundColor','white');
fprintf('Generated Figure S5:\n%s\n', outputFile)

end