sessionSD = movstd(allRewards, 10);
sessionSD = sessionSD(1:end-1);
[rwdSD_Sorted,tmpInds] = sort(sessionSD);
rwdSD_Inds = responseInds(tmpInds) + 1;

tercile = floor(length(rwdSD_Inds)/3);
rwdSDI_Inds = rwdSD_Inds(1:tercile);
rwdSDII_Inds = rwdSD_Inds(tercile+1:tercile*2);
rwdSDIII_Inds = rwdSD_Inds(tercile*2+1:end);


% rwdSDI_Inds = rwdSD_Inds(rwdSD_Inds == 0);
% rwdSDII_Inds = rwdSD_Inds(rwdSD_Inds == 0.447213595499958);
% rwdSDIII_Inds = rwdSD_Inds(rwdSD_Inds == 0.547722557505166);

rwdSDI_spikeMatx = [];
rwdSDII_spikeMatx = [];
rwdSDIII_spikeMatx = [];

for j = 1:length(rwdSDI_Inds)
    rwdSDI_spikeMatx(j,:) = allTrial_spikeMatx(rwdSDI_Inds(j),:);
end
for j = 1:length(rwdSDII_Inds)
    rwdSDII_spikeMatx(j,:) = allTrial_spikeMatx(rwdSDII_Inds(j),:);
end
for j = 1:length(rwdSDIII_Inds)
    rwdSDIII_spikeMatx(j,:) = allTrial_spikeMatx(rwdSDIII_Inds(j),:);
end


mySDF_rwdSDI = fastsmooth(nanmean(rwdSDI_spikeMatx, 1)*1000, 250);
mySDF_rwdSDII = fastsmooth(nanmean(rwdSDII_spikeMatx, 1)*1000, 250);
mySDF_rwdSDIII = fastsmooth(nanmean(rwdSDIII_spikeMatx, 1)*1000, 250);

figure; hold on;
plot(time, mySDF_rwdSDIII(1:length(time)),'b', 'LineWidth',2)
plot(time, mySDF_rwdSDII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
plot(time, mySDF_rwdSDI(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
legend({'High','Middle','Low'},'FontSize',5,'Location','northeast')