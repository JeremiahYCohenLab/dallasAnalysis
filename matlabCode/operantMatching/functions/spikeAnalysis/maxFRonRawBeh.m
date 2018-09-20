% time plot
rMag = 1;
nrMag = rMag/2;

figure; hold on
xlabel('Time (min)')
j = 1;
for i = 1:length(behSessionData)
    currTime = (behSessionData(i).CSon - behSessionData(1).CSon)/1000/60; %convert to min
    if strcmp(behSessionData(i).trialType,'CSplus')
        if ~isnan(behSessionData(i).rewardR)
            if behSessionData(i).rewardR == 1 % R side rewarded
                plot([currTime currTime],[0 rMag],'k')
            else
                plot([currTime currTime],[0 nrMag],'k') % R side not rewarded
            end
            plot([currTime currTime],[0 (maxFRtrial(j)/5)],'r')
            j=j+1;
        elseif ~isnan(behSessionData(i).rewardL)
            if behSessionData(i).rewardL == 1 % L side rewarded
                plot([currTime currTime],[-1*rMag 0],'k')
            else
                plot([currTime currTime],[-1*nrMag 0],'k')
            end
            plot([currTime currTime],[0 (maxFRtrial(j)/5)],'r')
            j=j+1;
            
        else % CSplus trial but no rewardL or rewardR
            plot([currTime currTime],[-rMag rMag],'r')
        end
    else % CS minus trial
        plot([currTime],0,'ko','markersize',4,'linewidth',2)
    end
    if any(i == origBlockSwitch)
        plot([currTime currTime],[-1*rMag rMag],'--','linewidth',1,'Color',[30 144 255]./255)
    end
end
xlim([0 currTime]);