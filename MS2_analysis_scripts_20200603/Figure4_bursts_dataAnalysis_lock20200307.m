%% Define colors
red = [230, 91, 98]/255;
BuOr = [33,102,172; 146,197,222; 253,174,97; 244,109,67] / 255;
YlPu = [255,242,135; 255,226,122; 254,195,105; 230,91,98; 205,98,126; 171,103,152; 134,106,180; 73,47,110] / 255;
BuPk = [38,63,182; 1,107,227; 235,70,87] / 255;
blue = [3,126,199]/255;

set(0, 'DefaultFigureColor', 'w', 'DefaultAxesFontSize', 7);

%% Load trajectories

load('M20171117_2min_dataMatrices_final.mat');

%%

% Substract minimum of p53 and MS2 for each single cell trajectory
MS2_normalization = min(traces_YFP')';
MS2_normalization = repmat(MS2_normalization, 1, size(traces_YFP,2));
traces_YFP = traces_YFP - MS2_normalization;

p53_normalization = min(traces_CFP')';
p53_normalization = repmat(p53_normalization, 1, size(traces_YFP,2));
traces_CFP = traces_CFP - p53_normalization;

% Plot the distribution of p21-MS2 intensity
figure; set(gcf, 'Position', [680 936 205 162]); [y,x] = hist(log(traces_YFP(traces_YFP > 0)), 100);
h = bar(x,y/sum(y)); set(h, 'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'none');
xlabel('log(p21-MS2 signal) [a.u.]'); ylabel('Relative frequency');
hold all; plot(([-0.48, -0.48]), ylim, 'r--', 'LineWidth', 1.5); xlim([-4,2]);

% Smooth p53 trajectories
w = 45;
blurred_CFP = smoothMatrix(traces_CFP, w, 'lowess');

peaks_CFP = getPeakMatrix_v3(blurred_CFP, 1, 5, 20);

%% Plot sample trajectories

set(gcf, 'Position', [560   793   560   155]);

selectedCell = find([annotation{:,2}] == 10 & [annotation{:,3}] == 49);
subplot(1,2,1); plotDualReporter((1:size(traces_YFP,2))/30 + 5, traces_CFP(selectedCell,:), traces_YFP(selectedCell,:), zeros(size(traces_YFP,2),1), [0, 100], [0,5], 'p53-CFP intensity (a.u.)', 'p21-MS2 signal (a.u.)', 'Time post-DNA damage (h)', [blue; [0.5, 0.5, 0.5]], gca);

%% Align trajectories around each p53 peak found previously

pastOffset = 300; futureOffset = 300; maxTimepoint = size(traces_CFP,2); timepoints = 1:maxTimepoint;

subCells = logical(ones(size(traces_YFP,1),1));
peakTiming = getDivisionTiming(peaks_CFP.traces_peaks_locs(:,timepoints));

aligned_CFP = alignTraces_TransitionEvents(blurred_CFP(subCells,timepoints), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_valleys_CFP = alignTraces_TransitionEvents(peaks_CFP.traces_valleys_locs(subCells,timepoints), peakTiming(subCells,:), pastOffset, futureOffset);
aligned_peaks_CFP = alignTraces_TransitionEvents(peaks_CFP.traces_peaks_locs(subCells,timepoints), peakTiming(subCells,:), pastOffset, futureOffset);

aligned_YFP = alignTraces_TransitionEvents(smoothMatrix(traces_YFP(subCells,timepoints), 5, 'movmedian'), peakTiming(subCells,:), pastOffset, futureOffset);

% Shift aligned p21-MS2 trajectories so that p21-MS2 and p53-CFP peaks
% coincide
shifting = 16;
aligned_YFP(:,shifting:end) = aligned_YFP(:,1:end-shifting+1); aligned_YFP(1:shifting-1)=NaN; 

% Find aligned peaks for which there is a right and left trough in the p53
% channel (complete pulses)
aligned_valid = NaN * ones(size(aligned_CFP));
for i=1:size(aligned_CFP,1)
    index1 = find(aligned_valleys_CFP(i,:) == 1 & (1:size(aligned_CFP,2)) < pastOffset-1, 1, 'last');
    if(~isempty(index1))
        aligned_valid(i,max(index1-1,1):end) = 1;
    end
    index2 = find(aligned_valleys_CFP(i,:) == 1 & (1:size(aligned_CFP,2)) > pastOffset+1, 1, 'first');
    if(~isempty(index2))
        aligned_valid(i,min(index2+1,size(aligned_CFP,2)):end) = NaN;
    end
end
aligned_valid(isnan(aligned_CFP)) = NaN;
validPeaks = nansum(aligned_valleys_CFP(:,1:(pastOffset-1)),2) == 1 & nansum(aligned_valleys_CFP(:,(pastOffset+1):end),2) == 1;

% Filter aligned trajectories to include only complete pulses
aligned_CFP = aligned_CFP(validPeaks,:);
aligned_valleys_CFP = aligned_valleys_CFP(validPeaks,:);
aligned_peaks_CFP = aligned_peaks_CFP(validPeaks,:);
aligned_YFP = aligned_YFP(validPeaks,:);
aligned_valid = aligned_valid(validPeaks,:);

%% Distribution of number of bursts per p53 peak

% Consider only the central peak in aligned matrices

ms2_threshold = exp(-0.376);
transitionMatrix = transitionEventMatrix(aligned_YFP, ms2_threshold, ms2_threshold, 1, 5, -1) .* aligned_valid;
numTransitions = nansum(transitionMatrix,2);
freq = tabulate(numTransitions);
figure; bar(freq(:,1), freq(:,3)/100); xlim([-1.05, 5.05]); ylim([0,1]);
xlabel('Number of bursts within a p53 pulse'); ylabel('Fraction of p53 pulses'); title('10 min');

%%
fraction_on = sum(aligned_YFP .* aligned_valid > ms2_threshold,2) ./ sum(~isnan(aligned_valid),2);

aligned_valid2 = NaN * ones(size(aligned_valid));
for i=1:size(aligned_valid,1)
    ms2_start = find(aligned_YFP(i,:) .* aligned_valid(i,:) > ms2_threshold, 1, 'first');
    ms2_end = find(aligned_YFP(i,:) .* aligned_valid(i,:) > ms2_threshold, 1, 'last');
    if(~isempty(ms2_start) && ~isempty(ms2_end))
        aligned_valid2(i,ms2_start:ms2_end) = 1;
    end
end

fraction_on2 = sum(aligned_YFP .* aligned_valid2 > ms2_threshold,2) ./ sum(~isnan(aligned_valid2),2);
fraction_on3 = fraction_on2;
fraction_on3(fraction_on3 == 1) = NaN;

figure;
boxplot(1- fraction_on3);
ylabel('% time p21-MS2 OFF');
ylim([0,1]);