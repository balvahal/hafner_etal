%% Read measurements table

red = [230, 91, 98]/255;
BuOr = [33,102,172; 146,197,222; 253,174,97; 244,109,67] / 255;
YlPu = [255,242,135; 255,226,122; 254,195,105; 230,91,98; 205,98,126; 171,103,152; 134,106,180; 73,47,110] / 255;
BuPk = [38,63,182; 1,107,227; 235,70,87] / 255;
blue = [3,126,199]/255;

set(0, 'DefaultAxesFontSize', 10, 'DefaultFigureColor', 'w');

%%
load('M20170724_biallelic_dataMatrices_final.mat');

%%

MS2_normalization = min(traces_YFP')';
MS2_normalization = repmat(MS2_normalization, 1, size(traces_YFP,2));
traces_YFP = traces_YFP - MS2_normalization;

blurred_YFP = smoothMatrix(traces_YFP, 3, 'movmedian');
blurred_CFP = smoothMatrix(traces_CFP, 8, 'movmean');

%%
ac_YFP = autocorrelationMatrix(smoothMatrix(blurred_YFP, 12, 'movmean'), 0:40);
ac_CFP = autocorrelationMatrix(smoothMatrix(blurred_CFP, 3, 'movmean'), 0:40);

peaks_ac_CFP = getPeakMatrix(ac_CFP, repmat(size(ac_CFP, 2), 1, size(ac_CFP,1))', 2);
peaks_ac_YFP = getPeakMatrix(ac_YFP, repmat(size(ac_YFP, 2), 1, size(ac_YFP,1))', 2);
pulsingCells = sum(peaks_ac_YFP.traces_peaks_locs(:,18:28),2) == 1 & sum(peaks_ac_YFP.traces_peaks_locs(:,1:15),2) == 0;

peaks_CFP = getPeakMatrix_v2(blurred_CFP, 1, 10);
peakTiming = getDivisionTiming(peaks_CFP.traces_peaks_locs);

minCFP = repmat(min(traces_CFP')', 1, size(traces_YFP,2));
normalized_CFP = traces_CFP - minCFP;
peakValues = normalized_CFP(peaks_CFP.traces_peaks_locs > 0);
normalizationFactor = mean(peakValues);
normalized_CFP = normalized_CFP ./ normalizationFactor;

normalized_Texas = traces_Texas;
maxTexas = exp(6.38); % Obtained as the mode of induced p21 in all cells: figure; hist(log(reshape(traces_Texas(:,1:96), [], 1)), 50)
normalized_Texas = normalized_Texas ./ maxTexas;

blurred_Texas = smoothMatrix(normalized_Texas, 8, 'movmean');
derivative_Texas = smoothMatrix(blurred_Texas(:,2:end) - blurred_Texas(:,1:end-1), 8, 'movmean');
derivative_Texas(:,end+1) = derivative_Texas(:,end);

%%
peakTiming_temp = peakTiming;
peakTiming_temp(peakTiming == 0) = NaN;

% Get unique cell ids. Each cell will have 2 alleles
freq = tabulate(cell_id(pulsingCells));
two_pulsing_cells = cellfun(@(x) x == 2, freq(:,2));
uniqueCellId = freq(two_pulsing_cells,1);

colocalization = zeros(3,length(uniqueCellId));
colocalization_p21 = zeros(3,length(uniqueCellId));

otherNuclei = zeros(2,length(uniqueCellId));

rng(1);

for i=1:length(uniqueCellId)
    % Find alleles of current cell
    currentCellIndex = find(strcmp(cell_id, uniqueCellId{i}));
    
    % Compute correlation between alleles in the same cell
    currentCells = blurred_YFP(currentCellIndex,:);
    colocalization(1,i) = corr(currentCells(1,:)', currentCells(2,:)', 'type', 'spearman');
    
    currentCells = derivative_Texas(currentCellIndex,:);
    colocalization_p21(1,i) = corr(currentCells(1,:)', currentCells(2,:)', 'type', 'spearman');
    
    % Identify alleles in randomly selected cell
    possible_cells = 1:length(cell_id);
    randomCellIndex = find(strcmp(cell_id, cell_id(randsample(possible_cells(~ismember(possible_cells,currentCellIndex)), 1))));
    
    % Compute correlations between an allele of the current cell and an
    % allele of a randomly selected cell
    currentCells = blurred_YFP(currentCellIndex,:);
    randomCells = blurred_YFP(randomCellIndex,:);
    colocalization(3,i) = corr(currentCells(1,:)', randomCells(2,:)', 'type', 'spearman');
    
    currentCells = derivative_Texas(currentCellIndex,:);
    randomCells = derivative_Texas(randomCellIndex,:);
    colocalization_p21(3,i) = corr(currentCells(1,:)', randomCells(2,:)', 'type', 'spearman');
    
    % Record the random nucleus that was chosen for the current analysis
    otherNuclei(2,i) = randomCellIndex(1);
    
    % Get peaks from current cell
    currentPeaks = peakTiming_temp(currentCellIndex(1),:);
    % Trim vector to consider only actual peaks
    currentPeaks = currentPeaks(~isnan(currentPeaks));
    % Generate matrix for comparison
    currentPeaks = repmat(currentPeaks', 1, size(peakTiming_temp,1))';
    % Get average distance between peaks for all other cells in the dataset
    peak_distance = sqrt((peakTiming_temp(:,1:size(currentPeaks,2)) - currentPeaks) .^ 2);
    peak_distance = sum(peak_distance,2);
    
    
    indexes = 1:length(peak_distance);
    minimumDistance = find(peak_distance == min(peak_distance(~ismember(indexes, currentCellIndex) & pulsingCells')));
    
    if(~isempty(minimumDistance))
        currentCells = blurred_YFP(currentCellIndex,:);
        randomCells = blurred_YFP(minimumDistance(1),:);
        colocalization(2,i) = corr(currentCells(1,:)', randomCells', 'type', 'spearman');
        
        currentCells = derivative_Texas(currentCellIndex,:);
        randomCells = derivative_Texas(minimumDistance(1),:);
        colocalization_p21(2,i) = corr(currentCells(1,:)', randomCells', 'type', 'spearman');
        otherNuclei(1,i) = minimumDistance(1);
    end
end

figure; set(gcf, 'Position', [680   560   285   200]);
boxplot(flipud(colocalization)', 'orientation', 'horizontal', 'Labels', {'unrelated random', 'unrelated in phase', 'biallelic'}); xlim([-1.05,1.05]);
xlabel('Spearman correlation'); title('p21-MS2');

ranksum(colocalization(1,:), colocalization(3,:))
ranksum(colocalization(2,:), colocalization(3,:))
ranksum(colocalization(1,:), colocalization(2,:))

%%
featured_cell = 'H2_IRCisNut1_53_4';

biallelic_cells_indexes = find(strcmp(cell_id, featured_cell));
currentCellIndex = find(strcmp(uniqueCellId, featured_cell));
inPhase_cells_indexes = otherNuclei(1,currentCellIndex(1));
outPhase_cell_indexes = otherNuclei(2,currentCellIndex(1));

%%
timepoints = (1:size(traces_YFP,2))/4;
figure; set(gcf, 'Position', [11   512   591   391]);
subplot(3,3,1); plot(timepoints, normalized_CFP(biallelic_cells_indexes,:)); xlim([0,24]); ylim([0,1.2]);
title('p53-CFP');
subplot(3,3,2); plot(timepoints, traces_YFP(biallelic_cells_indexes,:)); xlim([0,24]); ylim([0,3]);
title('p21-MS2')
subplot(3,3,4); plot(timepoints, normalized_CFP([biallelic_cells_indexes(1), inPhase_cells_indexes(1)],:)); xlim([0,24]); ylim([0,1.2]);
subplot(3,3,5); plot(timepoints, traces_YFP([biallelic_cells_indexes(1), inPhase_cells_indexes(1)],:)); xlim([0,24]); ylim([0,3]);
subplot(3,3,7); plot(timepoints, normalized_CFP([biallelic_cells_indexes(1), outPhase_cell_indexes],:)); xlim([0,24]); ylim([0,1.2]);
xlabel('Time post-DNA damage (h)');
subplot(3,3,8); plot(timepoints, traces_YFP([biallelic_cells_indexes(1), outPhase_cell_indexes],:)); xlim([0,24]); ylim([0,3]);
xlabel('Time post-DNA damage (h)');
