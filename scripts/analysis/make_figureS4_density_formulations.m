function make_figureS4_density_formulations(repoRoot)
%
% MAKE_FIGURES4_DENSITY_FORMULATIONS
%
% Generates Supporting Figure S4 comparing sea-ice density formulations
% used in CryoSat-2 sensitivity studies and this study.
%
% Formulations:
%   A2010 : Alexandrov et al. (2010), constant FYI/MYI densities
%   J2022 : Jutila et al. (2022), freeboard-dependent density
%   This study : freeboard- and snow-dependent limiting branches
%
% Output:
%   figures/FigS4.png
%

close all

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

%% Paths

figureDir = fullfile(repoRoot, ...
    'figures');

if ~exist(figureDir,'dir')
    mkdir(figureDir)
end

outputFigure = fullfile(figureDir, ...
    'FigS4.png');

%% Freeboard range

hf = linspace(0,0.4,500);   % sea-ice freeboard, m

%% Density formulations

% Alexandrov et al. (2010), constant density assumptions.
rho_A2010_FYI = 916.7;   % kg m^-3
rho_A2010_MYI = 882.0;   % kg m^-3

% Jutila et al. (2022), freeboard-dependent parameterization.
rho_J2022 = 72.0 .* exp(-3.7 .* hf) + 881.8;

% This study: low-snow / snow-free branch.
% Use the rounded continuous representation shown in Fig. 2.
hf_low_transition = 0.28;
rho_low = nan(size(hf));

idx_low = hf < hf_low_transition;
rho_low(idx_low)  = 865 + 125 .* hf(idx_low);
rho_low(~idx_low) = 900;

% This study: snow-covered branch, hs >= hscrit.
hf_high_transition = 0.20;
rho_snow = nan(size(hf));

idx_snow = hf <= hf_high_transition;
rho_snow(idx_snow)  = 911;
rho_snow(~idx_snow) = 911 - 20 .* (hf(~idx_snow) - hf_high_transition);

%% Plot settings

fs_ax = 10.5;
fs_txt = 9.0;
fsz_panel = 10;
lw = 2.0;

xlims = [0 0.4];
ylims = [860 960];

xticks_use = 0:0.1:0.4;
yticks_use = 860:20:960;

% Manuscript colors.
c_A2010 = [1 61 115] / 255;
c_this  = [58 174 140] / 255;
c_J2022 = [245 174 16] / 255;

%% Figure

fig = figure('Color','w');
set(fig,'Units','inches','Position',[1.1 1.0 10.5 3.15]);

tl = tiledlayout(fig,1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

ax = gobjects(1,3);

%% Panel a: A2010

ax(1) = nexttile(tl,1);
hold(ax(1),'on')

plot(ax(1),hf,rho_A2010_FYI .* ones(size(hf)),'-', ...
    'Color',c_A2010, ...
    'LineWidth',lw)

plot(ax(1),hf,rho_A2010_MYI .* ones(size(hf)),'--', ...
    'Color',c_A2010, ...
    'LineWidth',lw)

legend(ax(1),{'A2010 FYI','A2010 MYI'}, ...
    'Location','northeast', ...
    'Box','off', ...
    'FontSize',fs_txt, ...
    'Interpreter','tex')

title(ax(1),'A2010','FontWeight','normal','FontSize',fs_ax)
xlabel(ax(1),'Sea-ice freeboard {\it h_f} (m)','Interpreter','tex')
ylabel(ax(1),'Sea-ice density \rho_i (kg m^{-3})','Interpreter','tex')

xlim(ax(1),xlims)
ylim(ax(1),ylims)
xticks(ax(1),xticks_use)
yticks(ax(1),yticks_use)

set(ax(1), ...
    'FontSize',fs_ax, ...
    'FontWeight','normal', ...
    'LineWidth',1, ...
    'TickDir','in')

box(ax(1),'on')

%% Panel b: J2022

ax(2) = nexttile(tl,2);
hold(ax(2),'on')

plot(ax(2),hf,rho_J2022,'-', ...
    'Color',c_J2022, ...
    'LineWidth',lw)

legend(ax(2),{'J2022'}, ...
    'Location','northeast', ...
    'Box','off', ...
    'FontSize',fs_txt, ...
    'Interpreter','tex')

title(ax(2),'J2022','FontWeight','normal','FontSize',fs_ax)
xlabel(ax(2),'Sea-ice freeboard {\it h_f} (m)','Interpreter','tex')
ylabel(ax(2),'')

xlim(ax(2),xlims)
ylim(ax(2),ylims)
xticks(ax(2),xticks_use)
yticks(ax(2),yticks_use)

set(ax(2), ...
    'FontSize',fs_ax, ...
    'FontWeight','normal', ...
    'LineWidth',1, ...
    'TickDir','in')

box(ax(2),'on')

%% Panel c: This study

ax(3) = nexttile(tl,3);
hold(ax(3),'on')

plot(ax(3),hf,rho_snow,'-', ...
    'Color',c_this, ...
    'LineWidth',lw)

plot(ax(3),hf,rho_low,'--', ...
    'Color',c_this, ...
    'LineWidth',lw)

legend(ax(3),{'Snow-covered','Low snow'}, ...
    'Location','northeast', ...
    'Box','off', ...
    'FontSize',fs_txt, ...
    'Interpreter','tex')

title(ax(3),'This study','FontWeight','normal','FontSize',fs_ax)
xlabel(ax(3),'Sea-ice freeboard {\it h_f} (m)','Interpreter','tex')
ylabel(ax(3),'')

xlim(ax(3),xlims)
ylim(ax(3),ylims)
xticks(ax(3),xticks_use)
yticks(ax(3),yticks_use)

set(ax(3), ...
    'FontSize',fs_ax, ...
    'FontWeight','normal', ...
    'LineWidth',1, ...
    'TickDir','in')

box(ax(3),'on')

%% Panel labels

panelLabels = {'(a)','(b)','(c)'};
labelX = -0.16;
labelY = 1.05;

for i = 1:numel(ax)
    text(ax(i),labelX,labelY,panelLabels{i}, ...
        'Units','normalized', ...
        'FontSize',fsz_panel, ...
        'FontWeight','normal', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'Clipping','off')
end

set(findall(fig,'Type','axes'),'Toolbar',[])

%% Export

exportgraphics(fig,outputFigure, ...
    'Resolution',300)

close(fig)

fprintf('Generated Figure S4.\n')
fprintf('Saved figure to:\n%s\n',outputFigure)

end