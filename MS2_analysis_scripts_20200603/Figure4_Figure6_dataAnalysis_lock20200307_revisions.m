%% Define colors
red = [230, 91, 98]/255;
BuOr = [33,102,172; 146,197,222; 253,174,97; 244,109,67] / 255;
YlPu = [255,242,135; 255,226,122; 254,195,105; 230,91,98; 205,98,126; 171,103,152; 134,106,180; 73,47,110] / 255;
BuPk = [38,63,182; 1,107,227; 235,70,87] / 255;
blue = [3,126,199]/255;

set(0, 'DefaultFigureColor', 'w', 'DefaultAxesFontSize', 7);

%%
load('M20171129_dataMatrices_final.mat');

%%
p21_inducers = traces_Texas(:,60) > 20;
MS2_normalization = min(traces_YFP')';
MS2_normalization = repmat(MS2_normalization, 1, size(traces_YFP,2));
traces_YFP = traces_YFP - MS2_normalization;

p53_normalization = min(traces_CFP')';
p53_normalization = repmat(p53_normalization, 1, size(traces_CFP,2));
traces_CFP = traces_CFP - p53_normalization;

p21_normalization = min(traces_Texas')';
p21_normalization = repmat(p21_normalization, 1, size(traces_Texas,2));
traces_Texas = traces_Texas - p21_normalization;

blurred_YFP = smoothMatrix(traces_YFP, 3, 'movmedian');

[~, group_number] = ismember(annotation(:,1), unique(annotation(:,1)));

blurred_Texas = smoothMatrix(traces_Texas, 15, 'lowess');
derivative_Texas = smoothMatrix(blurred_Texas(:,2:end) - blurred_Texas(:,1:end-1), 15, 'lowess');
derivative_Texas(:,end+1) = derivative_Texas(:,end);

%%
timepoints = (1:size(traces_YFP,2))/12;

figure; set(gcf, 'Position', [178         278        1034         514]);
subplot(2,4,1); plotnfill_auto_quantiles_exclude(timepoints, (traces_CFP(group_number == 1,:)), 0.25, [0.3, 0.3, 0.3], -1); hold all; plotnfill_auto_quantiles_exclude(timepoints, (traces_CFP(group_number == 2,:)), 0.25, BuOr(1,:), -1);
xlabel('Time post irradiation (h)'); ylabel('p53-CFP intensity [a.u.]'); ylim([0,1000]);
subplot(2,4,2); plotnfill_auto_quantiles(timepoints, (traces_YFP(group_number == 1,:)), 0.25, [0.3, 0.3, 0.3]); hold all; plotnfill_auto_quantiles(timepoints, (traces_YFP(group_number == 2,:)), 0.25, BuOr(1,:));
xlabel('Time post irradiation (h)'); ylabel('p21-MS2 signal [a.u.]'); ylim([0,4]);
subplot(2,4,3); plotnfill_auto_quantiles(timepoints, (traces_Texas(group_number == 1 & p21_inducers,:)), 0.25, [0.3, 0.3, 0.3]); hold all; plotnfill_auto_quantiles(timepoints, (traces_Texas(group_number == 2 & p21_inducers,:)), 0.25, BuOr(1,:));
xlabel('Time post irradiation (h)'); ylabel('p21-mCherry intensity [a.u.]');
subplot(2,4,4); plotnfill_auto_quantiles_exclude(timepoints, (derivative_Texas(group_number == 1 & p21_inducers,:)), 0.25, [0.3, 0.3, 0.3], -1); hold all; plotnfill_auto_quantiles_exclude(timepoints, (derivative_Texas(group_number == 2 & p21_inducers,:)), 0.25, BuOr(1,:), -1);
xlabel('Time post irradiation (h)'); ylabel('p53-CFP intensity [a.u.]'); ylim([-2,10]);
subplot(2,4,5); plotnfill_auto_quantiles_exclude(timepoints, (traces_CFP(group_number == 1,:)), 0.25, [0.3, 0.3, 0.3], -1); hold all; plotnfill_auto_quantiles_exclude(timepoints, (traces_CFP(group_number == 2,:)), 0.25, BuOr(1,:), -1);
xlabel('Time post irradiation (h)'); ylabel('p53-CFP intensity [a.u.]'); xlim([0,10]); ylim([0,125]);
subplot(2,4,7); plotnfill_auto_quantiles(timepoints, (traces_Texas(group_number == 1 & p21_inducers,:)), 0.25, [0.3, 0.3, 0.3]); hold all; plotnfill_auto_quantiles(timepoints, (traces_Texas(group_number == 2 & p21_inducers,:)), 0.25, BuOr(1,:));
xlabel('Time post irradiation (h)'); ylabel('p21-mCherry intensity [a.u.]'); xlim([0,20]); ylim([0,100]);
subplot(2,4,8); plotnfill_auto_quantiles_exclude(timepoints(1:120), (derivative_Texas(group_number == 1 & p21_inducers,1:120)), 0.25, [0.3, 0.3, 0.3], -1); hold all; plotnfill_auto_quantiles_exclude(timepoints(1:120), (derivative_Texas(group_number == 2 & p21_inducers,1:120)), 0.25, BuOr(1,:), -1);
xlabel('Time post irradiation (h)'); ylabel('p53-CFP intensity [a.u.]'); ylim([-2,2]);

%% Smooth p53 trajectories
w=12;
blurred_CFP = smoothMatrix(traces_CFP, w, 'lowess');
peaks_CFP = getPeakMatrix_v3(blurred_CFP, 1, 10, 48);

derivative_CFP = smoothMatrix(smoothMatrix(blurred_CFP(:,2:end)+1, 1, 'movmean') - smoothMatrix(blurred_CFP(:,1:end-1)+1, 1, 'movmean'), 8, 'movmean');
derivative_CFP(:,end+1) = derivative_CFP(:,end);

%%
figure; set(gcf, 'Position', [680 936 205 162]); 
[y,x] = hist(reshape(log(blurred_YFP(group_number == 1,:)),[],1), 100);
h = bar(x,y/sum(y)); set(h, 'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'none');
xlabel('log(p21-MS2 signal) [a.u.]'); ylabel('Count'); xlim([-4,2]);
hold all; plot([-0.46, -0.46], ylim, 'r--', 'LineWidth', 1.5);

ms2_threshold = 0.63;

%% Align p53 trajectories on peaks
peakTiming = getDivisionTiming(peaks_CFP.traces_peaks_locs);
subCells = group_number == 1;
pastOffset = 40; futureOffset = 40;

id_matrix = repmat(1:size(traces_YFP,1)', size(traces_YFP,2), 1)';
timepoint_matrix = repmat(1:size(traces_YFP,2)', size(traces_YFP,1), 1);

aligned_CFP = alignTraces_TransitionEvents(blurred_CFP(subCells,:), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_CFP_derivative = alignTraces_TransitionEvents(derivative_CFP(subCells,:), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_YFP = alignTraces_TransitionEvents(blurred_YFP(subCells,:), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_YFP_blurred = alignTraces_TransitionEvents(smoothMatrix(blurred_YFP(subCells,:), 12, 'lowess'), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_id = alignTraces_TransitionEvents(id_matrix(subCells,:), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_timepoint = alignTraces_TransitionEvents(timepoint_matrix(subCells,:), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_derivative_Texas = alignTraces_TransitionEvents(derivative_Texas(subCells,:), peakTiming(subCells,:), pastOffset, futureOffset);

max_CFP = find(nanmean(aligned_CFP) == max(nanmean(aligned_CFP)));
max_YFP = find(nanmean(aligned_YFP_blurred) == max(nanmean(aligned_YFP_blurred)));

shifting = max_CFP - max_YFP + 1;

aligned_CFP_shifted = aligned_CFP(:,shifting:end);
aligned_CFP_derivative_shifted = aligned_CFP_derivative(:,shifting:end);
aligned_id_shifted = aligned_id(:,shifting:end);
aligned_timepoint_shifted = aligned_timepoint(:,shifting:end);
aligned_YFP_shifted = aligned_YFP(:,1:(end-shifting+1));
aligned_YFP_blurred_shifted = aligned_YFP_blurred(:,1:(end-shifting+1));
aligned_derivative_Texas_shifted = aligned_derivative_Texas(:,shifting:end);
blurred_CFP_shifted = blurred_CFP(:,shifting:end);
derivative_CFP_shifted = derivative_CFP(:,shifting:end);
blurred_YFP_shifted = blurred_YFP(:,1:(end-shifting+1));
derivative_Texas_shifted = derivative_Texas(:,shifting:end);
id_matrix_shifted = id_matrix(:,shifting:end);

%%

max_Texas = find(nanmean(aligned_derivative_Texas_shifted) == max(nanmean(aligned_derivative_Texas_shifted)));
max_CFP = find(nanmean(aligned_CFP_shifted) == max(nanmean(aligned_CFP_shifted)));
max_YFP = find(nanmean(aligned_YFP_blurred_shifted) == max(nanmean(aligned_YFP_blurred_shifted)));

shifting = max_Texas - max_YFP + 1;
aligned_CFP_shifted = aligned_CFP_shifted(:,1:(end-shifting+1));
aligned_CFP_derivative_shifted = aligned_CFP_derivative_shifted(:,1:(end-shifting+1));
aligned_id_shifted = aligned_id_shifted(:,1:(end-shifting+1));
aligned_timepoint_shifted = aligned_timepoint_shifted(:,1:(end-shifting+1));
aligned_YFP_shifted = aligned_YFP_shifted(:,1:(end-shifting+1));
aligned_YFP_blurred_shifted = aligned_YFP_blurred_shifted(:,1:(end-shifting+1));
blurred_CFP_shifted = blurred_CFP_shifted(:,1:(end-shifting+1));
derivative_CFP_shifted = derivative_CFP_shifted(:,1:(end-shifting+1));
blurred_YFP_shifted = blurred_YFP_shifted(:,1:(end-shifting+1));
aligned_derivative_Texas_shifted = aligned_derivative_Texas_shifted(:,shifting:end);
derivative_Texas_shifted = derivative_Texas_shifted(:,shifting:end);
id_matrix_shifted = id_matrix_shifted(:,shifting:end);

%%
% Input output relationship normalized

peak_loc = find(peaks_CFP.traces_peaks_locs);
timepoint_matrix = repmat(1:size(traces_YFP,2), size(traces_YFP,1), 1);

peak_value = peaks_CFP.traces_peaks_values(peak_loc);
peak_timepoint = timepoint_matrix(peak_loc);
peakNumber = discretizeValues(peak_timepoint, 22:22:96);
normalizationFactor = mean(aligned_CFP(:,pastOffset));

p53_binning = 0:0.05:2; xlimits = [-0.05,2];
discretized_p53 = discretizeValues(aligned_CFP_shifted(:) / normalizationFactor, p53_binning);
subValues = aligned_CFP_shifted(:) > 0 & aligned_YFP_shifted(:) > 0 & aligned_CFP_derivative_shifted(:) > 0;
mean_p53 = grpstats(aligned_CFP_shifted(subValues)  / normalizationFactor, discretized_p53(subValues), 'mean');
p_p21 = grpstats(aligned_YFP_shifted(subValues) > ms2_threshold, discretized_p53(subValues), 'mean');
numel = grpstats(aligned_YFP_shifted(subValues), discretized_p53(subValues), 'numel');
sem_p_p21 = grpstats(aligned_YFP_shifted(subValues) > ms2_threshold, discretized_p53(subValues), 'sem');
mean_p21 = grpstats(aligned_YFP_shifted(subValues & aligned_YFP_shifted(:) > ms2_threshold), discretized_p53(subValues & aligned_YFP_shifted(:) > ms2_threshold), 'mean');
sem_p21 = grpstats(aligned_YFP_shifted(subValues & aligned_YFP_shifted(:) > ms2_threshold), discretized_p53(subValues & aligned_YFP_shifted(:) > ms2_threshold), 'sem');
numel = grpstats(aligned_CFP_shifted(subValues), discretized_p53(subValues), 'numel');

minNum = 30;
valid = numel >= minNum;

figure; set(gcf, 'Position', [965   917   341   170]);
subplot(1,2,1); errorbar(mean_p53(valid), p_p21(valid), sem_p_p21(valid) * 1.96, 'o-', 'Color', [0.3, 0.3, 0.3]); hold all;
xlabel('p53-CFP intensity (a.u.)'); ylabel('Fraction active transcription'); ylim([0, 1.05]); xlim(xlimits);
subplot(1,2,2); errorbar(mean_p53(valid), mean_p21(valid), sem_p21(valid) * 1.96, 'o-', 'Color', [0.3, 0.3, 0.3]); hold all; 
xlabel('p53-CFP intensity (a.u.)'); ylabel('p21-MS2 signal (a.u.)'); ylim([0, 4]); xlim(xlimits);

%%
% Input output relationship normalized

peak_loc = find(peaks_CFP.traces_peaks_locs);
timepoint_matrix = repmat(1:size(traces_YFP,2), size(traces_YFP,1), 1);

peak_value = peaks_CFP.traces_peaks_values(peak_loc);
peak_timepoint = timepoint_matrix(peak_loc);
peakNumber = discretizeValues(peak_timepoint, 22:22:96);
normalizationFactor = mean(aligned_CFP(:,pastOffset));

p53_binning = 0:0.05:15; xlimits = [-0.05,15];
discretized_p53 = discretizeValues(aligned_CFP_shifted(:) / normalizationFactor, p53_binning);
subValues = aligned_CFP_shifted(:) > 0 & aligned_YFP_shifted(:) > 0 & aligned_CFP_derivative_shifted(:) > 0;
mean_p53 = grpstats(aligned_CFP_shifted(subValues)  / normalizationFactor, discretized_p53(subValues), 'mean');
p_p21 = grpstats(aligned_YFP_shifted(subValues) > ms2_threshold, discretized_p53(subValues), 'mean');
numel = grpstats(aligned_YFP_shifted(subValues), discretized_p53(subValues), 'numel');
sem_p_p21 = grpstats(aligned_YFP_shifted(subValues) > ms2_threshold, discretized_p53(subValues), 'sem');
mean_p21 = grpstats(aligned_YFP_shifted(subValues & aligned_YFP_shifted(:) > ms2_threshold), discretized_p53(subValues & aligned_YFP_shifted(:) > ms2_threshold), 'mean');
sem_p21 = grpstats(aligned_YFP_shifted(subValues & aligned_YFP_shifted(:) > ms2_threshold), discretized_p53(subValues & aligned_YFP_shifted(:) > ms2_threshold), 'sem');
numel = grpstats(aligned_CFP_shifted(subValues), discretized_p53(subValues), 'numel');

subValues1 = subValues;

minNum = 30;
valid = numel >= minNum;

figure; set(gcf, 'Position', [965   917   341   170]);
subplot(1,2,1); errorbar(mean_p53(valid), p_p21(valid), sem_p_p21(valid), 'o-', 'Color', [0.3, 0.3, 0.3]); hold all;
xlabel('p53-CFP intensity (a.u.)'); ylabel('Fraction active transcription'); ylim([0, 1.05]); xlim(xlimits);
subplot(1,2,2); errorbar(mean_p53(valid), mean_p21(valid), sem_p21(valid), 'o-', 'Color', [0.3, 0.3, 0.3]); hold all; 
xlabel('p53-CFP intensity (a.u.)'); ylabel('p21-MS2 signal (a.u.)'); ylim([0, 4]); xlim(xlimits);

%%

nutlin_CFP = blurred_CFP_shifted(group_number == 2,:);
nutlin_CFP_derivative = derivative_CFP_shifted(group_number == 2,:);
nutlin_YFP = blurred_YFP_shifted(group_number == 2,:);
annotation_nutlin = annotation(group_number == 2,:);

discretized_p53 = discretizeValues(nutlin_CFP(:) / normalizationFactor, p53_binning);
subValues = nutlin_CFP(:) > 0 & nutlin_YFP(:) > 0;
mean_p53 = grpstats(nutlin_CFP(subValues) / normalizationFactor, discretized_p53(subValues), 'mean');
p_p21 = grpstats(nutlin_YFP(subValues) > ms2_threshold, discretized_p53(subValues), 'mean');
mean_p21 = grpstats(nutlin_YFP(subValues & nutlin_YFP(:) > ms2_threshold), discretized_p53(subValues & nutlin_YFP(:) > ms2_threshold), 'mean');
sem_p21 = grpstats(nutlin_YFP(subValues & nutlin_YFP(:) > ms2_threshold), discretized_p53(subValues & nutlin_YFP(:) > ms2_threshold), 'sem');
numel = grpstats(nutlin_CFP(subValues), discretized_p53(subValues), 'numel');

valid = numel >= minNum;

subplot(1,2,1); errorbar(mean_p53(valid), p_p21(valid), sem_p21(valid), 'o-', 'Color', BuOr(1,:)); hold all; xlim(xlimits);
legend({'IR only', 'IR + MDM2i'}, 'location', 'southeast'); legend boxoff;
plot(repmat(nanmean(nutlin_CFP(subValues) / normalizationFactor),2,1), ylim, 'Color', BuOr(1,:));
plot(repmat(nanmean(aligned_CFP_shifted(subValues1) / normalizationFactor),2,1), ylim, 'Color', [0.3, 0.3, 0.3]);
subplot(1,2,2); errorbar(mean_p53(valid), mean_p21(valid), sem_p21(valid), 'o-', 'Color', BuOr(1,:)); hold all; xlim(xlimits);
plot(xlim, [ms2_threshold, ms2_threshold], '--', 'Color', [0.5, 0.5, 0.5]);
legend({'IR only', 'IR + MDM2i'}, 'location', 'southeast'); legend boxoff;
plot(repmat(nanmean(nutlin_CFP(subValues) / normalizationFactor),2,1), ylim, 'Color', BuOr(1,:));
plot(repmat(nanmean(aligned_CFP_shifted(subValues1) / normalizationFactor),2,1), ylim, 'Color', [0.3, 0.3, 0.3]);