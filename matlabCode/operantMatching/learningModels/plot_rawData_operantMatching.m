function plot_rawData_operantMatching(figHandle, sessionData, unCorrectedBlockSwitch, tMax)

subplot(figHandle); hold on
rMag = 1;
nrMag = rMag/2;

xlabel('Trials')
xlim([1 tMax])
ylim([-1 1])
for i = 1:length(sessionData)
    if strcmp(sessionData(i).trialType,'CSplus')
        if ~isnan(sessionData(i).rewardR)
            if sessionData(i).rewardR == 1 % R side rewarded
                plot([i i],[0 rMag],'k')
            else
                plot([i i],[0 nrMag],'k') % R side not rewarded
            end
        elseif ~isnan(sessionData(i).rewardL)
            if sessionData(i).rewardL == 1 % L side rewarded
                plot([i i],[-1*rMag 0],'k')
            else
                plot([i i],[-1*nrMag 0],'k')
            end

        else % CSplus trial but no rewardL or rewardR
            plot([i i],[-rMag rMag],'r')
        end
    else % CS minus trial
        plot(i,0,'ko','markersize',4,'linewidth',2)
    end
    if any(i == unCorrectedBlockSwitch)
        plot([i i],[-1*rMag rMag],'--','linewidth',1,'Color',[30 144 255]./255)
    end
end