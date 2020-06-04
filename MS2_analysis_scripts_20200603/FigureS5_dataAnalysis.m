load('STEL20200131_measurements.mat');

%%
MS2 = measurements.singleCellTracks_foci{2}(:,1) ./  measurements.singleCellTracks_median{2}(:,1);
MS2_median = measurements.singleCellTracks_median{2}(:,1);
MS2_foci = measurements.singleCellTracks_foci{2}(:,1);
FISH = measurements.singleCellTracks_foci{2}(:,2) ./  measurements.singleCellTracks_median{2}(:,2);
FISH_dilated = measurements.singleCellTracks_dilated{2}(:,2);
FISH_median = measurements.singleCellTracks_median{2}(:,2);
FISH_foci = measurements.singleCellTracks_foci{2}(:,2);
solidity = measurements.singleCellTracks_solidity(:,1);
area = measurements.singleCellTracks_area(:,1);
DAPI = measurements.singleCellTracks_integrated{1}(:,1);
annotation = measurements.cellAnnotation;
%%
validCells = solidity > 0.9;
MS2 = MS2(validCells);
MS2_foci = MS2_foci(validCells);
MS2_median = MS2_median(validCells);
FISH = FISH(validCells);
FISH_dilated = FISH_dilated(validCells);
FISH_foci = FISH_foci(validCells);
FISH_median = FISH_median(validCells);
annotation = annotation(validCells,:);

[~, group_number] = ismember(annotation(:,1), unique(annotation(:,1)));

%%
BuOr = [33,102,172; 146,197,222; 253,174,97; 244,109,67] / 255;

set(0, 'DefaultFigureColor', 'w')
figure; set(gcf, 'Position', [680   809   326   289]);
scatter(log(MS2(group_number == 1)), log(FISH(group_number == 1)), 20, BuOr(1,:)); box on; hold all;
scatter(log(MS2(group_number == 2)), log(FISH(group_number == 2)), 20, [0.5, 0.5, 0.5]); box on;
scatter(log(MS2([8,24])), log(FISH([8,24])), 20, 'r'); box on;
xlabel('log(p21-MS2 signal) [a.u.]'); ylabel('log(p21 FISH signal) [a.u.]');
xlim([0.25, 2.1]);
legend({'IR + MDM2i', 'IR only'}, 'location', 'northwest'); legend boxoff;

