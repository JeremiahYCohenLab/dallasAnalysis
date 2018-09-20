% Analysis_JC

tetID = 'T1';  
load([tetID '.mat']);   % each file contains time stapms for spikes (ts), number of trials (ntrials) and Odor onsets (OdorOn)

% Rater plot
figure;
subplot(5,1,1:2)
window = [-1 4];
trials = 1:ntrials;
raster12(ts,OdorOn,trials,window);      % plot raster
title(tetID)
axis([window 0+0.5 length(trials)+0.5]);
box off
ylabel('Trials')
% xlabel('Time - OdorOn (sec)')

% PSTH (Box filter)
tw = 20;  % time window for averaging (ms)
[y, t] = filt_PSTH(ts,OdorOn,trials,window,'Box',tw);
subplot(5,1,3)
m = mean(y);
e = std(y)/sqrt(size(y,1));
errorbar_patch(t,m,e,[.6 .6 .6]); hold on
plot(t,m,'k')
set(gca,'XLim',window)
box off
ax = axis; x = ax(1)+(ax(2)-ax(1))*0.6; y = ax(4)*0.4;
text(x,ax(4)*0.6,sprintf('Box filter (tau: %d ms)',tw))
% xlabel('Time - OdorOn (sec)')
ylabel('spikes/s')

% Gaussian filter
sigma = 10;  % standard deviation (ms)
[y, t] = filt_PSTH(ts,OdorOn,trials,window,'Gaussian',sigma);
subplot(5,1,4)
m = mean(y);
e = std(y)/sqrt(size(y,1));
errorbar_patch(t,m,e,[.6 .6 .6]); hold on
plot(t,m,'k')
set(gca,'XLim',window)
box off
ax = axis; x = ax(1)+(ax(2)-ax(1))*0.6; y = ax(4)*0.4;
text(x,ax(4)*0.6,sprintf('Gaussian filter (sigma: %d ms)',sigma))
% xlabel('Time - OdorOn (sec)')
ylabel('spikes/s')

% PSTH (EPSP filter)
tau = 20;   % time constant (ms)
[y, t] = filt_PSTH(ts,OdorOn,trials,window,'PSP',tau);
subplot(5,1,5)
m = mean(y);
e = std(y)/sqrt(size(y,1));
errorbar_patch(t,m,e,[.6 .6 .6]); hold on
plot(t,m,'k')
set(gca,'XLim',window)
box off
ax = axis; x = ax(1)+(ax(2)-ax(1))*0.6; y = ax(4)*0.4;
text(x,ax(4)*0.6,sprintf('PSP filter (tau: %d ms)',tau))
xlabel('Time - OdorOn (sec)')
ylabel('spikes/s')
