 figure; hold on
    mySDF_rwdR = fastsmooth(nanmean(rwd_spikeMatx(lickRrwd_Inds,:), 1)*1000, smoothWin);
    mySDF_rwdL = fastsmooth(nanmean(rwd_spikeMatx(lickLrwd_Inds,:), 1)*1000, smoothWin);
    mySDF_noRwdR = fastsmooth(nanmean(noRwd_spikeMatx(lickRnoRwd_Inds,:), 1)*1000, smoothWin);
    mySDF_noRwdL = fastsmooth(nanmean(noRwd_spikeMatx(lickLnoRwd_Inds,:), 1)*1000, smoothWin);
    plot(time, mySDF_rwdR(1:length(time)),'-','Color', [0.5 0 0.8],'LineWidth',2)
    plot(time, mySDF_noRwdR(1:length(time)),'-','Color',[0.9 0.6 1],'LineWidth',2)
    plot(time, mySDF_rwdL(1:length(time)),'-b','LineWidth',2)
    plot(time, mySDF_noRwdL(1:length(time)),'-','Color', [0.6 0.6 1],'LineWidth',2)
%    legend({'rwd R','no rwd R','rwd L','no rwd L'},'FontSize',14,'Location','northeast')
    ylabel('Spikes/s')
    xlabel('Time (ms)')
    xlim([-500 2000])
    set(gca, 'linewidth',2,'fontsize', 14)
    
    
    
    
figure; hold on;
    mySDF_rwdPrev = fastsmooth(nanmean(rwdPrev_spikeMatx, 1)*1000, smoothWin);
    mySDF_noRwdPrev = fastsmooth(nanmean(noRwdPrev_spikeMatx, 1)*1000, smoothWin);
    plot(time, mySDF_rwdPrev(1:length(time)),'-b','LineWidth',2)
    plot(time, mySDF_noRwdPrev(1:length(time)),'-','Color',[0.7 0 1],'LineWidth',2)
    legend('rwd', 'no rwd')
    ylabel('Spikes/s')
    xlabel('Time (ms)')
    xlim([-500 2000])
    set(gca, 'linewidth',2,'fontsize', 14)
    
    
 smoothWin = 250; 
figure; 
    subplot(1,2,1); hold on
    mySDF_rwdHxI = fastsmooth(nanmean(rwdHxI_spikeMatx, 1)*1000, smoothWin);
    mySDF_rwdHxII = fastsmooth(nanmean(rwdHxII_spikeMatx, 1)*1000, smoothWin);
    mySDF_rwdHxIII = fastsmooth(nanmean(rwdHxIII_spikeMatx, 1)*1000, smoothWin);
    plot(time, mySDF_rwdHxIII(1:length(time)),'b', 'LineWidth',2)
    plot(time, mySDF_rwdHxII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
    plot(time, mySDF_rwdHxI(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
    legend({'High','Middle','Low'},'FontSize',5,'Location','northeast')
    ylabel('Spikes/s')
    xlabel('Time (ms)')
    xlim([-500 2000])
    set(gca, 'linewidth',2,'fontsize', 14)
    
    subplot(1,2,2); hold on;
    plotSpikeRaster(allTrial_spike(i,rwdHx_Inds),'PlotType','vertline'); hold on
    set(gca,'Xticklabel',[],'linewidth',2,'fontsize', 14);
    xlabel('Time (ms)')
    ylabel('Trials')
    