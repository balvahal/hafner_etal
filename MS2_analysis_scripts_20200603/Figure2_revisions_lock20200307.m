%% Read measurements table

red = [230, 91, 98]/255;
BuOr = [33,102,172; 146,197,222; 253,174,97; 244,109,67] / 255;
YlPu = [255,242,135; 255,226,122; 254,195,105; 230,91,98; 205,98,126; 171,103,152; 134,106,180; 73,47,110] / 255;
BuPk = [38,63,182; 1,107,227; 235,70,87] / 255;
blue = [3,126,199]/255;

set(0, 'DefaultAxesFontSize', 7, 'DefaultFigureColor', 'w');

%%
load('M20170724_dataMatrices_final.mat');

%%

% Substract the minimum value from each p53 and p21-MS2 trajectory, preserving the
% identity of missing datapoints
w=8;
blurred_CFP = smoothMatrix(traces_CFP, w, 'lowess');
for i=1:size(traces_CFP,1)
    traces_CFP(i,:) = traces_CFP(i,:) - min(blurred_CFP(i,traces_CFP(i,:) > 0));
    traces_CFP(i,traces_CFP(i,:) < 0) = 0;
    
    traces_YFP(i,:) = traces_YFP(i,:) - min(traces_YFP(i,traces_YFP(i,:) > 0));
    traces_YFP(i,traces_YFP(i,:) < 0) = 0;
end


%%
% Obtain the timing of p21 induction and p21 degradation from each cell

p21_threshold = exp(4.7);
p21_dropTiming = zeros(size(traces_Texas,1),1);
p21_inductionTiming = zeros(size(traces_Texas,1),1);
validMatrix = zeros(size(traces_Texas));
for i=1:length(p21_inductionTiming)
    currentTrace = traces_Texas(i,:) > p21_threshold;
    inductionTiming = findpattern_once(currentTrace, ones(5,1));
    if(~isempty(inductionTiming))
        validMatrix(i,inductionTiming:end) = 1;        
        p21_inductionTiming(i) = inductionTiming;
        dropTiming = findpattern_once(currentTrace(inductionTiming:end), zeros(5,1));
        if(~isempty(dropTiming))
            p21_dropTiming(i) = dropTiming + inductionTiming - 1;
            validMatrix(i,p21_dropTiming(i):end) = 0;
        end
    end
end

% Identify cells which induce p21 within the first p53 pulse (inducer
% cells) and cells which start with high p21 values (initiallyInduced)
validCells = validMatrix(:,20) == 1 & validMatrix(:,end) == 1;
inducerCells = validMatrix(:,1) == 0 & validCells;
initiallyInduced = validMatrix(:,1) == 1 & validCells;

%%

% Smooth single cell trajectories
w=8;
blurred_YFP = smoothMatrix(traces_YFP, w, 'lowess');
blurred_CFP = smoothMatrix(traces_CFP, w, 'lowess');
blurred_Texas = smoothMatrix(traces_Texas, w, 'lowess');
baseValueMatrix = repmat(blurred_CFP(:,1), 1, size(blurred_CFP,2));

% Estimate p21 protein derivative
derivative_Texas = smoothMatrix(smoothMatrix(blurred_Texas(:,2:end)+1, 8, 'lowess') - smoothMatrix(blurred_Texas(:,1:end-1)+1, 8, 'lowess'), 8, 'lowess');
derivative_Texas(:,end+1) = derivative_Texas(:,end);

% Call peaks in p53, p21-MS2 and p21 protein derivative
peaks_CFP = getPeakMatrix_v3(blurred_CFP, 1, 40, 10);
peaks_derivative_Texas = getPeakMatrix_v3(derivative_Texas, 1, 10, 16);
peaks_YFP = getPeakMatrix_v3(blurred_YFP, 1, 0.25, 10);

% Plot summary statistics of the different molecular species
figure; set(gcf, 'Position', [62   883   795   187]);
timepoints = (1:size(traces_Texas,2))/4;
subplot(1,4,1); plotnfill_auto_quantiles(timepoints, blurred_CFP(validCells,:), 0.25, blue); ylim([0,750]); xlabel('Time post-irradiation (h)'); ylabel('p53-CFP intensity (a.u.)');
subplot(1,4,2); plotnfill_auto_quantiles(timepoints, blurred_YFP(validCells,:), 0.25, [0.4, 0.4, 0.4]); ylim([0,2.5]); xlabel('Time post-irradiation (h)'); ylabel('p21-MS2 signal (a.u.)');
subplot(1,4,3); plotnfill_auto_quantiles(timepoints, blurred_Texas(validCells,:), 0.25, BuOr(4,:)); xlabel('Time post-irradiation (h)'); ylabel('p21-mCherry intensity (a.u.)');
subplot(1,4,4); plotnfill_auto_quantiles(timepoints, derivative_Texas(validCells,:), 0.25, BuOr(4,:)); xlabel('Time post-irradiation (h)'); ylabel('p21-mCherry derivative (a.u.)');

[~, group_number] = ismember(annotation(:,1), unique(annotation(:,1)));

%%

% Calculate cross-correlations between all pairs of molecular species
xcorr_p53_ms2 = crosscorrelationMatrix(blurred_YFP(validCells,:), blurred_CFP(validCells,:), 0:40);
xcorr_ms2_p21 = crosscorrelationMatrix(blurred_YFP(validCells,:), derivative_Texas(validCells,:), 0:40);
xcorr_p53_p21 = crosscorrelationMatrix(blurred_CFP(validCells,:), derivative_Texas(validCells,:), 0:40);

% Calculate autocorrelations of p53, p21-MS2 and p21-mCherry derivative
xcorr_p53_p53 = crosscorrelationMatrix(blurred_CFP(validCells,:), blurred_CFP(validCells,:), 0:40);
xcorr_ms2_ms2 = crosscorrelationMatrix(blurred_YFP(validCells,:), blurred_YFP(validCells,:), 0:40);
xcorr_p21_p21 = crosscorrelationMatrix(derivative_Texas(validCells,:), derivative_Texas(validCells,:), 0:40);

% Plot cross-correlations and autocorrelations
timepoints = (0:(size(xcorr_ms2_p21,2)-1))/4;
figure; set(gcf, 'Position', [680 953 560 145]);
subplot(1,3,1); 
plotnfill_auto_quantiles(timepoints, xcorr_p53_p53, 0.25, blue); ylim([-0.5,1]); hold all;
plotnfill_auto_quantiles(timepoints, xcorr_ms2_ms2, 0.25, [0.5, 0.5, 0.5]); ylim([-0.5,1]); hold all;
plotnfill_auto_quantiles(timepoints, xcorr_p21_p21, 0.25, BuOr(4,:)); ylim([-0.5,1]); hold all;
ylabel('Autocorrelation'); xlabel('Time lag (h)');

%%

% Define parameters for aligning single cell trajectories on the basis of
% p53 pulses
pastOffset = 20; futureOffset = 20; maxTimepoint = 96; timepoints = 1:maxTimepoint;
id_matrix = repmat(1:size(traces_YFP,1), size(traces_YFP,2), 1)';

% Consider all cells for alignment. You can modify subCells to consider
% only a subpopulation of cells
subCells = logical(ones(size(traces_YFP,1),1));
% Obtain timing of p53 peaks in all trajectories
peakTiming = getDivisionTiming(peaks_CFP.traces_peaks_locs(:,timepoints));

% Align trajectories on the basis of p53 peak timing
aligned_CFP = alignTraces_TransitionEvents(blurred_CFP(subCells,timepoints), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_valleys_CFP = alignTraces_TransitionEvents(peaks_CFP.traces_valleys_locs(subCells,timepoints), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_peaks_CFP = alignTraces_TransitionEvents(peaks_CFP.traces_peaks_locs(subCells,timepoints), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_YFP = alignTraces_TransitionEvents(traces_YFP(subCells,timepoints), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_YFP_blurred = alignTraces_TransitionEvents(blurred_YFP(subCells,timepoints), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_Texas = alignTraces_TransitionEvents(blurred_Texas(subCells,timepoints), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_derivative_Texas = alignTraces_TransitionEvents(derivative_Texas(subCells,timepoints), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_id = alignTraces_TransitionEvents(id_matrix(subCells,timepoints), peakTiming(subCells,:), pastOffset, futureOffset);

% Identify peaks for which there are left and right valleys (full pulses)
validPeaks = nansum(aligned_valleys_CFP(:,1:(pastOffset-1)),2) == 1 & nansum(aligned_valleys_CFP(:,(pastOffset+1):end),2) == 1;

% Include only full pulses for further analyses
aligned_derivative_Texas = aligned_derivative_Texas(validPeaks,:);
aligned_CFP = aligned_CFP(validPeaks,:);
aligned_valleys_CFP = aligned_valleys_CFP(validPeaks,:);
aligned_peaks_CFP = aligned_peaks_CFP(validPeaks,:);
aligned_YFP = aligned_YFP(validPeaks,:);
aligned_YFP_blurred = aligned_YFP_blurred(validPeaks,:);
aligned_id = aligned_id(validPeaks,:);

aligned_valid = NaN * ones(size(aligned_CFP));
for i=1:size(aligned_CFP,1)
    indexes = find(aligned_valleys_CFP(i,:) == 1);
    aligned_valid(i,indexes(1):indexes(2)) = 1;
end

%% 
% Shift aligned trajectories to account for time delays in peak p21-MS2 and
% p21 protein derivative in relationship to p53 pulses.

timepoints = (-pastOffset:futureOffset) / 4;

figure; set(gcf, 'Position', [560   517   206   431]);
subplot(3,1,1); plot(timepoints, nanmean(aligned_CFP), 'LineWidth', 2); 
ylabel('Average p53-CFP intensity (a.u.)'); xlabel('Time from p53 peak (h)');
subplot(3,1,2); plot(timepoints, nanmean(aligned_YFP), 'LineWidth', 2); hold all;
ylabel('Average p21-MS2 intensity (a.u.)'); xlabel('Time from p53 peak (h)');

shifting = 3;
aligned_YFP(:,shifting:end) = aligned_YFP(:,1:end-shifting+1); aligned_YFP(1:shifting-1)=NaN; 
aligned_YFP_blurred(:,shifting:end) = aligned_YFP_blurred(:,1:end-shifting+1); aligned_YFP_blurred(1:shifting-1)=NaN; 
subplot(3,1,2); plot(timepoints, nanmean(aligned_YFP), 'LineWidth', 2); legend({'Before shifting', 'After shifting'}); legend boxoff;
subplot(3,1,3); plot(timepoints, nanmean(aligned_derivative_Texas), 'LineWidth', 2); hold all;

shifting = 5;
aligned_Texas(:,1:end-shifting+1) = aligned_Texas(:,shifting:end); aligned_Texas(:,(end-shifting+1):end)=NaN; 
aligned_derivative_Texas(:,1:end-shifting+1) = aligned_derivative_Texas(:,shifting:end); aligned_derivative_Texas(:,(end-shifting+1):end)=NaN; 
subplot(3,1,3); plot(timepoints, nanmean(aligned_derivative_Texas), 'LineWidth', 2); legend({'Before shifting', 'After shifting'}); legend boxoff;
ylabel('Average p21-mCherry derivative (a.u.)'); xlabel('Time from p53 peak (h)');

%%

subCells = ismember(aligned_id(:,pastOffset), find(inducerCells));
figure; set(gcf, 'Position', [560   789   368   300]);
subplot(2,2,1);
xval = reshape(aligned_YFP_blurred(subCells,:),1,[]);
yval = reshape(aligned_derivative_Texas(subCells,:), 1, []);
validPoints = ~isnan(xval) & ~isnan(yval) & xval > 0;
scatter((xval(validPoints)'), yval(validPoints)', '.');
ylim([-150,200]); box on; xlabel('p21-MS2 intensity (a.u.)'); ylabel('p21-mCherry derivative (a.u.)');
text(2, 200, sprintf('Spearman r = %0.2f', corr(xval(validPoints)', yval(validPoints)', 'type', 'spearman')))

subplot(2,2,2);
xval = reshape(aligned_YFP_blurred(subCells,:),1,[]);
yval = reshape(aligned_Texas(subCells,:), 1, []);
validPoints = ~isnan(xval) & ~isnan(yval);
scatter(xval(validPoints)', yval(validPoints)', '.');
box on; xlabel('p21-MS2 intensity (a.u.)'); ylabel('p21-mCherry intensity (a.u.)');
text(2, 6000, sprintf('Spearman r = %0.2f', corr(xval(validPoints)', yval(validPoints)', 'type', 'spearman')))

subplot(2,2,3);
xval = log(aligned_CFP(:));
yval = log(aligned_YFP(:));
validPoints = ~isnan(xval) & ~isnan(yval);
scatter(xval(validPoints)', yval(validPoints)', '.');
ylim([-6,3]); box on; xlabel('log(p53-CFP intensity) [a.u.]'); ylabel('log(p21-MS2 intensity) [a.u.]');
text(5, 2, sprintf('Spearman r = %0.2f', corr(xval(validPoints), yval(validPoints), 'type', 'spearman')))

