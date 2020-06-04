load('M20171118_dataMatrix_fixVsLive_final.mat');

set(0, 'DefaultFigureColor', 'w');

%%
validCells = traces_solidity(:,end) > 0.9 & traces_Texas(:,end) > exp(3);
validCells = traces_Texas(:,end) > exp(3);

figure; set(gcf, 'Position', [560 610 555 338]);
subplot(2,3,1); scatterCorrelation(log(traces_CFP(validCells,end)), log(FISH(validCells))); xlim([2, 7]); xlabel('log(p53-CFP) [a.u.]'); ylim([4.5, 8]); ylabel('log(p21 smFISH) [a.u.]'); box on;
subplot(2,3,2); scatterCorrelation(log(FISH(validCells,end)), log(traces_Texas(validCells,end))); xlabel('log(p21 smFISH) [a.u.]'); ylim([1.7250, 6.5804]); ylabel('log(p21-mCherry) [a.u.]'); box on;
subplot(2,3,3); scatterCorrelation(log(traces_CFP(validCells,end)), log(traces_Texas(validCells,end))); xlim([2, 7]); ylim([1.7250, 6.5804]); xlabel('log(p53-CFP) [a.u.]'); ylabel('log(p21-mCherry) [a.u.]'); box on;

%%
figure;
set(gcf, 'Position', [560 805 166 143])
scatterCorrelation(log(FISH(validCells,end)), log(traces_YFP(validCells,end)));
xlabel('ln(FISH signal) [a.u.]');
ylabel('ln(MS2 signal) [a.u.]');
xlim([3.5, 8]); ylim([0,2.5]);
box on;