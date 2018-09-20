figure
plotSpikeRaster(allTrial_spike(i,rBar_Inds),'PlotType','vertline'); hold on
set(gca,'Xticklabel',[]);
set(gca, 'xtick',[0 500]);
set(gca,'Yticklabel',[]);
set(gca,'linewidth',6)
set(gca, 'FontSize', 18)


half = floor(length(rBar_Inds)/2);
rBarI_Inds = rBar_Inds(1:half);
rBarII_Inds = rBar_Inds(half+1:end);


rBarI_spikeMatx = NaN(length(rBarI_Inds), length(time));
rBarII_spikeMatx = NaN(length(rBarII_Inds), length(time));

for j = 1:length(rBarI_Inds)
    rBarI_spikeMatx(j,:) = allTrial_spikeMatx(rBarI_Inds(j),:);
end
for j = 1:length(rBarII_Inds)
    rBarII_spikeMatx(j,:) = allTrial_spikeMatx(rBarII_Inds(j),:);
end

mySDF_rBarI = fastsmooth(nanmean(rBarI_spikeMatx, 1)*1000, 400);
mySDF_rBarII = fastsmooth(nanmean(rBarII_spikeMatx, 1)*1000, 400);
mySDF_rBarIII = fastsmooth(nanmean(rBarIII_spikeMatx, 1)*1000, 400);

figure; hold on;
plot(time, mySDF_rBarI(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
plot(time, mySDF_rBarII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
plot(time, mySDF_rBarIII(1:length(time)),'b', 'LineWidth',2)

trialBeg = tb*1000;
CSoff = tb*1000 + 500;
for j = 1:length(allTrial_spike)
    if ~isempty(allTrial_spikeMatx(i,j))
        postCSspikeCount(j) = sum(allTrial_spikeMatx(j, trialBeg:CSoff));
    else
        postCSspikeCount(i,j) = NaN;
    end
end

