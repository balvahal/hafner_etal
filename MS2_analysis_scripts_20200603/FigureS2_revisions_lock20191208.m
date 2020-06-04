%% Read measurements table

red = [230, 91, 98]/255;
BuOr = [33,102,172; 146,197,222; 253,174,97; 244,109,67] / 255;
YlPu = [255,242,135; 255,226,122; 254,195,105; 230,91,98; 205,98,126; 171,103,152; 134,106,180; 73,47,110] / 255;
BuPk = [38,63,182; 1,107,227; 235,70,87] / 255;
blue = [3,126,199]/255;

set(0, 'DefaultAxesFontSize', 7, 'DefaultFigureColor', 'w');

%%
load('M20170724_dataMatrices_multipleDoses.mat');

%%

% Smooth single cell trajectories
w=8;
blurred_CFP = smoothMatrix(traces_CFP, w, 'lowess');

[~, group_number] = ismember(annotation(:,1), unique(annotation(:,1)));

%%
blues_gradient = flipud(cbrewer('seq', 'Blues', 9));
timepoints = (1:size(traces_CFP,2))/4;

figure; set(gcf, 'Position', [62   883   200   187]);
selectedGroup = 1;
additionalFilter = group_number == selectedGroup;
plotnfill_auto_quantiles(timepoints, blurred_CFP(additionalFilter,:), 0.25, blues_gradient(4,:)); ylim([0,750]); xlabel('Time post-irradiation (h)'); ylabel('p53-CFP intensity (a.u.)'); hold all;

%figure; set(gcf, 'Position', [62   883   795   187]);
selectedGroup = 2;
additionalFilter = group_number == selectedGroup;
plotnfill_auto_quantiles(timepoints, blurred_CFP(additionalFilter,:), 0.25, blues_gradient(6,:)); ylim([0,750]); xlabel('Time post-irradiation (h)'); ylabel('p53-CFP intensity (a.u.)'); hold all;

% Plot summary statistics of the different molecular species
selectedGroup = 3;
additionalFilter = group_number == selectedGroup;
plotnfill_auto_quantiles(timepoints, blurred_CFP(additionalFilter,:), 0.25, blues_gradient(8,:)); ylim([0,750]); xlabel('Time post-irradiation (h)'); ylabel('p53-CFP intensity (a.u.)'); hold all;