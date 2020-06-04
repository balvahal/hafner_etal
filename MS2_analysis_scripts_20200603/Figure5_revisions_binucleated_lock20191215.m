%% Read measurements table

red = [230, 91, 98]/255;
BuOr = [33,102,172; 146,197,222; 253,174,97; 244,109,67] / 255;
YlPu = [255,242,135; 255,226,122; 254,195,105; 230,91,98; 205,98,126; 171,103,152; 134,106,180; 73,47,110] / 255;
BuPk = [38,63,182; 1,107,227; 235,70,87] / 255;
blue = [3,126,199]/255;

set(0, 'DefaultAxesFontSize', 10, 'DefaultFigureColor', 'w');

%%
load('M20171211_binucleated_dataMatrix_final.mat');

%%
MS2_normalization = min(traces_YFP')';
MS2_normalization = repmat(MS2_normalization, 1, size(traces_YFP,2));
traces_YFP = traces_YFP - MS2_normalization;

% Use a median filter of 15 min to remove noise from p21-MS2 trajectories
blurred_YFP = smoothMatrix(traces_YFP, 3, 'movmean');
blurred_CFP = smoothMatrix(traces_CFP, 24, 'lowess');

%% Obtain pulsing cells

% Compute autocorrelation functions for each p53 and p21-MS2 trajectory
ac_YFP = autocorrelationMatrix(smoothMatrix(blurred_YFP, 12, 'movmean'), 0:160);
ac_CFP = autocorrelationMatrix(smoothMatrix(blurred_CFP, 3, 'movmean'), 0:160);

% Identify peaks in p53 autocorrelations and get pulsing cells as those
% which have a peak at a time delay between 60 and 100 frames, and a valley 
% between 1 and 50 frames
peaks_ac_CFP = getPeakMatrix(ac_CFP, repmat(size(ac_CFP, 2), 1, size(ac_CFP,1))', 2);
pulsingCells = sum(peaks_ac_CFP.traces_peaks_locs(:,60:100),2) == 1 & sum(peaks_ac_CFP.traces_peaks_locs(:,1:50),2) == 0;

%% Obtain the peaks of p53 trajectories

peaks_CFP = getPeakMatrix_v2(blurred_CFP, 2, 10);
peakTiming = getDivisionTiming(peaks_CFP.traces_peaks_locs);

% Normalize p53 trajectories. Substract the minimum value from each
% trajectory and divide by the mean peak value
minCFP = repmat(min(traces_CFP')', 1, size(traces_YFP,2));
normalized_CFP = traces_CFP - minCFP;
normalizationFactor = repmat(max(normalized_CFP')', 1, size(traces_YFP,2));
peakValues = normalized_CFP(peaks_CFP.traces_peaks_locs > 0);
normalizationFactor = mean(peakValues);
normalized_CFP = normalized_CFP ./ normalizationFactor;

% Normalize p21 trajectories. Substract the minimum value from each
% trajectory and divide by the the mode of the p21 intensity distribution
normalized_Texas = traces_Texas;
maxTexas = exp(5.2); % Obtained from the distribution of induced p21 in all cells
normalized_Texas = normalized_Texas ./ maxTexas;

% Smooth p21 protein trajectories and estimate derivatives
blurred_Texas = smoothMatrix(normalized_Texas, 24, 'movmean');
derivative_Texas = smoothMatrix(blurred_Texas(:,2:end) - blurred_Texas(:,1:end-1), 24, 'movmean');
derivative_Texas(:,end+1) = derivative_Texas(:,end);

%% Identify binucleated and mononucleated cells

% Cell id is the cell from which trajectories came from
cell_id = strcat(annotation(:,1), '_', cellfun(@num2str, annotation(:,2), 'UniformOutput', 0));
cell_freq = tabulate(cell_id(pulsingCells,:));
% Binucleated cells are those in which a cell id is repeated exactly 2
% times. Mononucleated cells are those in which a cell id is unique
binucleated_cells = cell_freq([cell_freq{:,2}] == 2,1);
mononucleated_cells = cell_freq([cell_freq{:,2}] == 1,1);

% Obtain the indexes of cells belonging to either the mononucleated or
% binucleated classes
binucleated_cell_idx = ismember(cell_id, binucleated_cells);
mononucleated_cell_idx = ismember(cell_id, mononucleated_cells);

%% p53 in binucleated cell analysis

cell_freq = tabulate(cell_id);
binucleated_cells_all = cell_freq([cell_freq{:,2}] == 2,1);

% What is the correlation between pairs of p53 trajectories?
p53_correlation = NaN * ones(length(binucleated_cells_all),2);

for i=1:length(p53_correlation)
    currentCells = find(strcmp(cell_id, binucleated_cells_all{i}));
    p53_correlation(i,1) = corr(traces_CFP(currentCells(1),:)', traces_CFP(currentCells(2),:)');
    
    randomCell = randsample(find(~strcmp(cell_id, binucleated_cells_all{i})), 1);
    p53_correlation(i,2) = corr(traces_CFP(currentCells(1),:)', traces_CFP(randomCell,:)');
end

figure; set(gcf, 'Position', [560   767   164   181]);
boxplot(p53_correlation, 'Labels', {'binucleated', 'unrelated'}); ylabel('p53 trajectory correlation'); ylim([-1.05,1.05]);

%% Plot examples of correlated and uncorrelated cells

% We selected the following cells to show the distinct ways in which p53
% trajectories are uncorrelated
uncorrelatedCells = {'M20171211_20_8', 'M20171211_20_5', 'M20171211_20_22', 'M20171211_3_8'};
timepoints = (1:size(traces_YFP,2))/12;

figure; set(gcf, 'Position', [680   530   453   568]);
for i=1:length(uncorrelatedCells)
    binucleated_cells_indexes = strcmp(cell_id, uncorrelatedCells{i});
    subplot(length(uncorrelatedCells),2,(i-1)*2+1); plot(timepoints, smoothMatrix(traces_CFP(binucleated_cells_indexes,:), 7, 'movmedian')); xlim([0,24]);
    subplot(length(uncorrelatedCells),2,(i-1)*2+2); plot(timepoints, traces_YFP(binucleated_cells_indexes,:)); xlim([0,24]); ylim([0,4]);
end

rng(80)
binucleatedPulsing = ismember(binucleated_cells_all, binucleated_cells);

correlatedCells = binucleated_cells_all(randsample(find(p53_correlation(:,1) > 0.6 & binucleatedPulsing), length(uncorrelatedCells)));
figure; set(gcf, 'Position', [680   530   453   568]);
for i=1:length(uncorrelatedCells)
    binucleated_cells_indexes = strcmp(cell_id, correlatedCells{i});
    subplot(length(uncorrelatedCells),2,(i-1)*2+1); plot(timepoints, smoothMatrix(traces_CFP(binucleated_cells_indexes,:), 7, 'movmedian')); xlim([0,24]);
    subplot(length(uncorrelatedCells),2,(i-1)*2+2); plot(timepoints, traces_YFP(binucleated_cells_indexes,:)); xlim([0,24]); ylim([0,4]);
end

%% Get the correlation of p21-MS2 and p21 protein trajectories in pairs of binucleated or unrelated cells
peakTiming_temp = peakTiming;
peakTiming_temp(peakTiming == 0) = NaN;

colocalization = zeros(3,length(binucleated_cells));
colocalization_p21 = zeros(3,length(binucleated_cells));

rng(1);

for i=1:length(binucleated_cells)
    % Identify which cells belong to the current binucleated pair
    currentCellIndex = find(strcmp(cell_id, binucleated_cells{i}));
    
    % Compute correlation between binucleated trajectories
    currentCells = blurred_YFP(currentCellIndex,:);
    colocalization(1,i) = corr(currentCells(1,:)', currentCells(2,:)', 'type', 'Spearman');
    colocalization_p21(1,i) = corr(derivative_Texas(currentCellIndex(1),:)', derivative_Texas(currentCellIndex(2),:)', 'type', 'Spearman');
    
    % Sample a random cell from the rest of the dataset, excluding both
    % nuclei from binucleated cell
    possible_cells = 1:length(binucleated_cells);
    randomCellIndex = find(strcmp(cell_id, binucleated_cells(randsample(possible_cells(possible_cells ~= i), 1))));

    % Compute correlation between unrelated cells
    randomCells = blurred_YFP(randomCellIndex,:);
    colocalization(3,i) = corr(currentCells(1,:)', randomCells(2,:)', 'type', 'Spearman');   
    colocalization_p21(3,i) = corr(derivative_Texas(currentCellIndex(1),:)', derivative_Texas(randomCellIndex(1),:)', 'type', 'Spearman');
    
    % Identify the peaks in p53 dynamics in the binucleated cells. Get an
    % unrelated cell with similar p53 profile, excluding both nuclei from
    % binucleated cell
    currentPeaks = peakTiming_temp(currentCellIndex(1),:);
    currentPeaks = currentPeaks(~isnan(currentPeaks));
    currentPeaks = repmat(currentPeaks', 1, size(peakTiming_temp,1))';
    peak_distance = sqrt((peakTiming_temp(:,1:size(currentPeaks,2)) - currentPeaks) .^ 2);
    peak_distance = sum(peak_distance,2);
    
    indexes = 1:length(peak_distance);
    minimumDistance = find(peak_distance == min(peak_distance(~ismember(indexes, currentCellIndex) & pulsingCells')));
    
    if(~isempty(minimumDistance))
        % Calculate correlations between pairs of nuclei of unrelated cells
        % which nonetheless share similar p53 trajectories.
        randomCells = blurred_YFP(minimumDistance(1),:);
        colocalization(2,i) = corr(currentCells(1,:)', randomCells', 'type', 'Spearman');
        colocalization_p21(2,i) = corr(derivative_Texas(currentCellIndex(1),:)', derivative_Texas(minimumDistance(1),:)', 'type', 'Spearman');
    end
end

figure; set(gcf, 'Position', [680   560   285   200]);
boxplot(flipud(colocalization)', 'orientation', 'horizontal', 'Labels', {'unrelated random', 'unrelated in phase', 'binucleated'}); xlim([-1.05,1.05]);
xlabel('Spearman correlation'); title('p21-MS2');

ranksum(colocalization(1,:), colocalization(3,:))
ranksum(colocalization(2,:), colocalization(3,:))
ranksum(colocalization(1,:), colocalization(2,:))

%% Plot sample pairs of binucleated and unrelated cells (with similar p53 trajectories or randomly sampled)

binucleated_cells_indexes = find(strcmp(cell_id, 'M20171211_11_5'));
inPhase_cells_indexes = find(strcmp(cell_id, 'M20171211_13_9'));
outPhase_cell_indexes = 18;

timepoints = (1:size(traces_YFP,2))/12;
figure; set(gcf, 'Position', [11   512   591   374]);
subplot(3,3,1); plot(timepoints, normalized_CFP(binucleated_cells_indexes,:)); xlim([0,24]); ylim([0,1.2]);
subplot(3,3,2); plot(timepoints, blurred_YFP(binucleated_cells_indexes,:)); xlim([0,24]); ylim([0,3]);

subplot(3,3,4); plot(timepoints, normalized_CFP([binucleated_cells_indexes(1), inPhase_cells_indexes(1)],:)); xlim([0,24]); ylim([0,1.2]);
subplot(3,3,5); plot(timepoints, blurred_YFP([binucleated_cells_indexes(1), inPhase_cells_indexes(1)],:)); xlim([0,24]); ylim([0,3]);

subplot(3,3,7); plot(timepoints, normalized_CFP([binucleated_cells_indexes(1), outPhase_cell_indexes],:)); xlim([0,24]); ylim([0,1.2]);
xlabel('Time post-DNA damage (h)');
subplot(3,3,8); plot(timepoints, blurred_YFP([binucleated_cells_indexes(1), outPhase_cell_indexes],:)); xlim([0,24]); ylim([0,3]);
xlabel('Time post-DNA damage (h)');