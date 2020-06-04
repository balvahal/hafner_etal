%%
red = [230, 91, 98]/255;
BuOr = [33,102,172; 146,197,222; 253,174,97; 244,109,67] / 255;
YlPu = [255,242,135; 255,226,122; 254,195,105; 230,91,98; 205,98,126; 171,103,152; 134,106,180; 73,47,110] / 255;
%spectral = vertcat(flipdim(cbrewer('seq','GnBu',100),1), cbrewer('seq', 'YlOrRd', 100));
BuPk = [38,63,182; 1,107,227; 235,70,87] / 255;
blue = [3,126,199]/255;

set(0, 'DefaultFigureColor', 'w', 'DefaultAxesFontSize', 12);

%%
load('M20181214_cellCycle_dataMatrix_final.mat');

%%
divisionTiming = getDivisionTiming(divisions);

[~, ordering] = sort(divisionTiming(:,1));

damageAdded = 97;

traces_CFP = traces_CFP(ordering,:);
traces_YFP = traces_YFP(ordering,:);
traces_Texas = traces_Texas(ordering,:);
traces_area = traces_area(ordering,:);
annotation = annotation(ordering,:);
divisions = divisions(ordering,:);
divisionTiming = divisionTiming(ordering,:);

for i=1:size(traces_YFP,1)
    traces_YFP(i,:) = traces_YFP(i,:) - min(traces_YFP(i,:));
    traces_CFP(i,:) = traces_CFP(i,:) - min(traces_CFP(i,:));
    traces_Texas(i,:) = traces_Texas(i,:) - min(traces_Texas(i,:));
    traces_CFP(i,:) = traces_CFP(i,:) / max(traces_CFP(i,1:damageAdded));
end

%%
[y,x] = hist(log(reshape(traces_CFP(:,1:damageAdded), 1, [])), 100);
geminin_threshold = exp(-1.73);
figure; set(gcf, 'Position', [560   834   155   114]);
set(gca, 'FontSize', 7);
bar(x,y/sum(y)); xlim([-8,0]); xlabel('ln(CFP-hGeminin(1-100) intensity) [a.u.]'); ylabel('Relative frequency');
hold all; plot(log([geminin_threshold, geminin_threshold]), ylim);
[y,x] = hist(log(reshape(traces_Texas(:,damageAdded:end), 1, [])), 100);
p21_threshold = exp(0.302);
figure; set(gcf, 'Position', [560   834   155   114]);
set(gca, 'FontSize', 7);
bar(x,y/sum(y)); xlim([-5,5]); xlabel('ln(p21-mCherry intensity) [a.u.]'); ylabel('Relative frequency');
hold all; plot(log([p21_threshold, p21_threshold]), ylim);

responded = traces_Texas(:,damageAdded + 5*4) > p21_threshold;
geminin_positive = traces_CFP(:,damageAdded) > geminin_threshold;

g1_cells = responded & ~geminin_positive;
g2_cells = responded & geminin_positive;
s_cells = ~responded & geminin_positive;
g1_s_cells = ~responded & ~geminin_positive;

cell_cycle_class = zeros(size(traces_CFP,1),1);
cell_cycle_class(g1_cells) = 0;
cell_cycle_class(g2_cells) = 3;
cell_cycle_class(s_cells) = 2;
cell_cycle_class(g1_s_cells) = 1;

[~, ordering] = sort(cell_cycle_class);
traces_CFP = traces_CFP(ordering,:);
traces_YFP = traces_YFP(ordering,:);
traces_Texas = traces_Texas(ordering,:);
traces_area = traces_area(ordering,:);
annotation = annotation(ordering,:);
divisions = divisions(ordering,:);
divisionTiming = divisionTiming(ordering,:);
cell_cycle_class = cell_cycle_class(ordering);

%%
figure; set(gcf, 'Position', [555   556   777   230]);
timepoints = [-fliplr(1:96), 0:(size(traces_YFP,2)-97)]/96;
subplot(1,3,1); plotnfill_auto_quantiles(timepoints, traces_CFP, 0.25, blue); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title('CFP-geminin'); ylabel('CFP-hGeminin(1-100) intensity [a.u.]')
subplot(1,3,2); plotnfill_auto_quantiles(timepoints, traces_YFP, 0.25, [0.4, 0.4, 0.4]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title('p21-MS2'); ylabel('p21-MS2 signal [a.u.]')
subplot(1,3,3); plotnfill_auto_quantiles(timepoints, traces_Texas, 0.25, BuOr(4,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title('p21-mKate2'); ylabel('p21-mKate2 intensity [a.u.]')

figure; set(gcf, 'Position', [555   556   777   230]);
timepoints = [-fliplr(1:96), 0:(size(traces_YFP,2)-97)]/96;
subCells = cell_cycle_class == 0; keyword = ' G1 cells';
subplot(1,3,1); plotnfill_auto_quantiles(timepoints, traces_CFP(subCells,:), 0.25, blue); ylim([0,1.5]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('CFP-geminin', keyword)); ylabel('CFP-hGeminin(1-100) intensity [a.u.]')
subplot(1,3,2); plotnfill_auto_quantiles(timepoints, traces_YFP(subCells,:), 0.25, [0.4, 0.4, 0.4]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('p21-MS2', keyword)); ylabel('p21-MS2 signal [a.u.]')
subplot(1,3,3); plotnfill_auto_quantiles(timepoints, traces_Texas(subCells,:), 0.25, BuOr(4,:));ylim([0,20]);   hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('p21-mCherry', keyword)); ylabel('p21-mKate2 intensity [a.u.]')

figure; set(gcf, 'Position', [555   556   777   230]);
timepoints = [-fliplr(1:96), 0:(size(traces_YFP,2)-97)]/96;
subCells = cell_cycle_class == 1; keyword = ' G1-S cells';
subplot(1,3,1); plotnfill_auto_quantiles(timepoints, traces_CFP(subCells,:), 0.25, blue); ylim([0,1.5]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('CFP-geminin', keyword)); ylabel('CFP-hGeminin(1-100) intensity [a.u.]')
subplot(1,3,2); plotnfill_auto_quantiles(timepoints, traces_YFP(subCells,:), 0.25, [0.4, 0.4, 0.4]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('p21-MS2', keyword)); ylabel('p21-MS2 signal [a.u.]')
subplot(1,3,3); plotnfill_auto_quantiles(timepoints, traces_Texas(subCells,:), 0.25, BuOr(4,:));ylim([0,20]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('p21-mCherry', keyword)); ylabel('p21-mKate2 intensity [a.u.]')

figure; set(gcf, 'Position', [555   556   777   230]);
timepoints = [-fliplr(1:96), 0:(size(traces_YFP,2)-97)]/96;
subCells = cell_cycle_class == 2; keyword = ' S cells';
subplot(1,3,1); plotnfill_auto_quantiles(timepoints, traces_CFP(subCells,:), 0.25, blue); ylim([0,1.5]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('CFP-geminin', keyword)); ylabel('CFP-hGeminin(1-100) intensity [a.u.]')
subplot(1,3,2); plotnfill_auto_quantiles(timepoints, traces_YFP(subCells,:), 0.25, [0.4, 0.4, 0.4]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('p21-MS2', keyword)); ylabel('p21-MS2 signal [a.u.]')
subplot(1,3,3); plotnfill_auto_quantiles(timepoints, traces_Texas(subCells,:), 0.25, BuOr(4,:)); ylim([0,20]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('p21-mCherry', keyword)); ylabel('p21-mKate2 intensity [a.u.]')

figure; set(gcf, 'Position', [555   556   777   230]);
timepoints = [-fliplr(1:96), 0:(size(traces_YFP,2)-97)]/96;
subCells = cell_cycle_class == 3; keyword = ' G2 cells';
subplot(1,3,1); plotnfill_auto_quantiles(timepoints, traces_CFP(subCells,:), 0.25, blue); ylim([0,1.5]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('CFP-geminin', keyword)); ylabel('CFP-hGeminin(1-100) intensity [a.u.]')
subplot(1,3,2); plotnfill_auto_quantiles(timepoints, traces_YFP(subCells,:), 0.25, [0.4, 0.4, 0.4]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('p21-MS2', keyword)); ylabel('p21-MS2 signal [a.u.]')
subplot(1,3,3); plotnfill_auto_quantiles(timepoints, traces_Texas(subCells,:), 0.25, BuOr(4,:)); ylim([0,20]); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(strcat('p21-mCherry', keyword)); ylabel('p21-mKate2 intensity [a.u.]')

figure; set(gcf, 'Position', [555   556   777   230]);
timepoints = [-fliplr(1:96), 0:(size(traces_YFP,2)-97)]/96;
subCells = cell_cycle_class == 3;
subplot(1,3,1); plotnfill_auto_quantiles(timepoints, traces_CFP(subCells,:), 0.5, YlPu(8,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('CFP-geminin')); ylabel('CFP-hGeminin(1-100) intensity [a.u.]'); hold all;
subplot(1,3,2); plotnfill_auto_quantiles(timepoints, traces_YFP(subCells,:), 0.5, YlPu(8,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('p21-MS2')); ylabel('p21-MS2 signal [a.u.]'); hold all;
subplot(1,3,3); plotnfill_auto_quantiles(timepoints, traces_Texas(subCells,:), 0.5, YlPu(8,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('p21-mCherry')); ylabel('p21-mKate2 intensity [a.u.]'); hold all;
subCells = cell_cycle_class == 0;
subplot(1,3,1); plotnfill_auto_quantiles(timepoints, traces_CFP(subCells,:), 0.5, YlPu(3,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('CFP-geminin')); ylabel('CFP-hGeminin(1-100) intensity [a.u.]'); hold all;
subplot(1,3,2); plotnfill_auto_quantiles(timepoints, traces_YFP(subCells,:), 0.5, YlPu(3,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('p21-MS2')); ylabel('p21-MS2 signal [a.u.]'); hold all;
subplot(1,3,3); plotnfill_auto_quantiles(timepoints, traces_Texas(subCells,:), 0.5, YlPu(3,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('p21-mCherry')); ylabel('p21-mKate2 intensity [a.u.]'); hold all;
subCells = cell_cycle_class == 1;
subplot(1,3,1); plotnfill_auto_quantiles(timepoints, traces_CFP(subCells,:), 0.5, YlPu(5,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('CFP-geminin')); ylabel('CFP-hGeminin(1-100) intensity [a.u.]'); hold all;
subplot(1,3,2); plotnfill_auto_quantiles(timepoints, traces_YFP(subCells,:), 0.5, YlPu(5,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('p21-MS2')); ylabel('p21-MS2 signal [a.u.]'); hold all;
subplot(1,3,3); plotnfill_auto_quantiles(timepoints, traces_Texas(subCells,:), 0.5, YlPu(5,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('p21-mCherry')); ylabel('p21-mKate2 intensity [a.u.]'); hold all;
subCells = cell_cycle_class == 2;
subplot(1,3,1); plotnfill_auto_quantiles(timepoints, traces_CFP(subCells,:), 0.5, YlPu(6,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('CFP-geminin')); ylabel('CFP-hGeminin(1-100) intensity [a.u.]'); hold all;
subplot(1,3,2); plotnfill_auto_quantiles(timepoints, traces_YFP(subCells,:), 0.5, YlPu(6,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('p21-MS2')); ylabel('p21-MS2 signal [a.u.]'); hold all;
subplot(1,3,3); plotnfill_auto_quantiles(timepoints, traces_Texas(subCells,:), 0.5, YlPu(6,:)); hold all; plot([0, 0], ylim, 'k--', 'LineWidth', 2); xlabel('Time post-irradiation (days)'); title(('p21-mCherry')); ylabel('p21-mKate2 intensity [a.u.]'); hold all;

g1_cells = find(cell_cycle_class == 0);
s_cells = find(cell_cycle_class == 2);
g2_cells = find(cell_cycle_class == 3);