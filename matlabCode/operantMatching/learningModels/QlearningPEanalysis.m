filename = [sessionName '.asc'];

[pe_Sorted,peInds] = sort(pe);
%posLine = find(outputStruct.pe > 0, 1);

figure; hold on
suptitle(sessionName)

%sort rasters by prediction error
subplot(2,4,1)
figure;
if ~all(cellfun(@isempty,allTrial_spike(i,responseInds(peInds))))
    plotSpikeRaster(allTrial_spike(i,responseInds(peInds)),'PlotType','vertline'); hold on
end
title('prediction error: low to high')


%sort rasters by local prediction error variance
peSD = movstd(pe, [5 0]);
[peSD_Sorted,peSDinds] = sort(peSD);

subplot(2,4,2)
if ~all(cellfun(@isempty,allTrial_spike(i,responseInds(peSDinds))))
    plotSpikeRaster(allTrial_spike(i,responseInds(peSDinds)),'PlotType','vertline'); hold on
end
title('prediction error variance: low to high')



%sort rasters by deviation from local average
for j = 5:length(pe)
    tmpZ = zscore(pe(j-4:j));
    peZ(j-4) = tmpZ(end);
end
[peZ_Sorted,peZinds] = sort(peZ);
peZinds = peZinds + 4;

subplot(2,4,3)
if ~all(cellfun(@isempty,allTrial_spike(i,responseInds(peZinds))))
    plotSpikeRaster(allTrial_spike(i,responseInds(peZinds)),'PlotType','vertline'); hold on
end
title('prediction error error: low to high')




%get pre-CS spike counts
    for j = 1:size(response_spikeMatx,1)
        preCSspikeCounts(j) = sum(response_spikeMatx(j, 1:trialBeg));
    end





for j = 1:length(peInds)
    peInds_spikeMatx(j,:) = allTrial_spikeMatx(responseInds(peSDinds(j)),:);
end

mySDF_peI = fastsmooth(nanmean(peInds_spikeMatx(1:100,:), 1)*1000, 250);
mySDF_peII = fastsmooth(nanmean(peInds_spikeMatx(101:200,:), 1)*1000, 250);
mySDF_peIII = fastsmooth(nanmean(peInds_spikeMatx(201:300,:), 1)*1000, 250);
figure; hold on
plot(time, mySDF_peI(1:length(time)),'-b','LineWidth',2)
plot(time, mySDF_peII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
plot(time, mySDF_peIII(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)

for j = 1:length(pe)
    trialTimes(j) = behSessionData(responseInds(j)).CSon - behSessionData(1).CSon;
end



figure; hold on; 
yyaxis right; plot(sessionSpikeSDF)
plot(sessionSpikeSDFsmooth,'-b')
ylim([2 7])
yyaxis left; plot(trialTimes, pe)
ylim([-0.1 0.6])



figure; hold on;
subplot(2,1,1)
yyaxis left; plot(peSD); 
ylim([0 1]); 
ylabel('Q learning prediction error variance (6 trials)')
yyaxis right; plot(preCSspikeCounts); 
ylabel('pre CS spike counts')
title('mCG16d20170529 TT4 SS02')
subplot(2,1,2)
crosscorr(preCSspikeCounts, peSD)


%linear regression analysis
    for j = 1:size(response_spikeMatx,1)
        CSspikeCounts(j) = sum(response_spikeMatx(j, trialBeg:(trialBeg+500)));
    end

tbl = table(CSspikeCounts(2:end)',rwdHx(1:end-1)',pe(1:end-1)',peSD(1:end-1)','VariableNames',{'CSspikeCounts','rwdHx','pe','peSD'});
meow = fitlm(tbl,'CSspikeCounts~rwdHx+pe+peSD')

%% rwd rate

rwdTimes = [sessionData(responseInds(2:end)).CSon] - sessionData(1).CSon;
rwdRate = sessionRwdsSmooth(rwdTimes);
figure; yyaxis left; plot(preCSspikeCounts); yyaxis right; plot(rwdRate)

[rwdRate_Sorted, rwdRateInds] = sort(rwdRate);

figure;
if ~all(cellfun(@isempty,allTrial_spike(i,responseInds(rwdRateInds))))
    plotSpikeRaster(allTrial_spike(i,responseInds(rwdRateInds)),'PlotType','vertline'); hold on
end

for j = 1:length(rwdRateInds)
    rwdRateInds_spikeMatx(j,:) = allTrial_spikeMatx(responseInds(rwdRateInds(j)),:);
end

mySDF_rwdRateI = fastsmooth(nanmean(rwdRateInds_spikeMatx(1:100,:), 1)*1000, 250);
mySDF_rwdRateII = fastsmooth(nanmean(rwdRateInds_spikeMatx(101:200,:), 1)*1000, 250);
mySDF_rwdRateIII = fastsmooth(nanmean(rwdRateInds_spikeMatx(201:300,:), 1)*1000, 250);
figure; hold on
plot(time, mySDF_rwdRateI(1:length(time)),'-b','LineWidth',2)
plot(time, mySDF_rwdRateII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
plot(time, mySDF_rwdRateIII(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)