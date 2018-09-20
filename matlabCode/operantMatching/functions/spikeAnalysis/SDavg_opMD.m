session_spikeArray2 = session_spikeArray(1:length(session_rwdsArray));
sessionSD = movstd(session_rwdsArray, 30000);
sessionMean = movmean(session_rwdsArray, 30000);
sessionSpike = movmean(session_spikeArray2, 30000);
figure; hold on;

purp = [0.5 0 0.8];
blue = [0 0 1];
set(gcf,'defaultAxesColorOrder',[blue; purp]);

subplot(2,2,1)
yyaxis right; plot(sessionMean,'-b'); 
ylabel('smoothed rewards')
yyaxis left; plot(sessionSpike, '-', 'Color', purp)
ylabel('smoothed spikes')
xlabel('time (ms)')

subplot(2,2,3)
scatter(sessionMean, sessionSpike);


subplot(2,2,2)
yyaxis left; plot(sessionSD, '-', 'Color', purp); 
ylabel('standard deviation of rewards')
yyaxis right; plot(sessionSpike, '-b')
ylabel('smoothed spikes')
xlabel('time (ms)')


subplot(2,2,4)
scatter(sessionSD, sessionSpike);