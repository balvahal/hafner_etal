function [] = plotDualReporter(time, YFP, RFP, division, ylim1, ylim2, ylab1, ylab2, xlab, color, axisHandle)
    red = [237,28,36] / 255;
    green = [0.1, 0.7, 0.1];
    blue = [3,126,199]/255;
    yellow = [255,203,5] / 255;
    %yellow = [255,181,5] / 255;
    %yellow = green;
    
    fontSize = 12;
    
    if(isempty(ylim1))
        ylim1 = [min(YFP) * 0.95, max(YFP) * 1.05];
    end
    if(isempty(ylim2))
        ylim2 = [min(RFP) * 0.95, max(RFP) * 1.05];
    end
    
    %f = figure; set(gcf, 'Position', [50,50,600,450], 'Units', 'Pixels', 'Color', 'w');
    %[ax,p1,p2] = plotyy(time,YFP,time,RFP,@stairs,@stairs);
    [ax,p1,p2] = plotyy(axisHandle, time,YFP,time,RFP);
    set(p1, 'Color', color(1,:), 'LineWidth', 1.5);
    set(p2, 'Color', color(2,:), 'LineWidth', 1.5);
    set(ax(1), 'FontSize', fontSize,'ycolor', color(1,:)); set(ax(2), 'FontSize', fontSize, 'ycolor', color(2,:));
    xlabel(ax(1), xlab, 'FontSize', fontSize); 
    ylabel(ax(1), ylab1, 'FontSize', fontSize); 
    ylabel(ax(2), ylab2, 'FontSize', fontSize);
    set(ax(1), 'ylim', ylim1, 'YTick', floor([min(ylim1):(diff(ylim1)/5):max(ylim1)]));
    set(ax(2), 'ylim', ylim2, 'YTick', ([min(ylim2):(diff(ylim2)/5):max(ylim2)]));
    set(ax(1), 'ylim', ylim1, 'xlim', [min(time), max(time)]);
    set(ax(2), 'ylim', ylim2, 'xlim', [min(time), max(time)]);
    
    hold all;
    divisionEvents = find(division > 0);
    for i =1:length(divisionEvents)
        plot(ax(1), time([divisionEvents(i), divisionEvents(i)]), ylim1, 'Color', [0.3, 0.3, 0.3]);
    end
end