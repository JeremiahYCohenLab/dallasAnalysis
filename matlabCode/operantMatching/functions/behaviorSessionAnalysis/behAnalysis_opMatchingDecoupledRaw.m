function [allRewards, allChoices] = behAnalysis_opMatchingDecoupledRaw(filename, saveFigFlag)

if nargin < 2
    saveFigFlag = 1;
end


[root, sep] = currComputer();

[animalName, date] = strtok(filename, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];
sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep filename '_sessionData_behav.mat'];
if exist(sessionDataPath,'file')
    load(sessionDataPath)
else
    [sessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(filename);
end

%% Break session down into CS+ trials where animal responded

responseInds = find(~isnan([sessionData.rewardTime])); % find CS+ trials with a response in the lick window
omitInds = isnan([sessionData.rewardTime]); 

tempBlockSwitch = blockSwitch;
for i = 2:length(blockSwitch)
    subVal = sum(omitInds(tempBlockSwitch(i-1):tempBlockSwitch(i)));
    blockSwitch(i:end) = blockSwitch(i:end) - subVal;
end

allReward_R = [sessionData(responseInds).rewardR]; 
allReward_L = [sessionData(responseInds).rewardL]; 
rewProb_R = [sessionData(responseInds).rewardProbR]; 
rewProb_L = [sessionData(responseInds).rewardProbL]; 
allChoices = NaN(1,length(sessionData(responseInds)));
allChoices(~isnan(allReward_R)) = 1;
allChoices(~isnan(allReward_L)) = -1;

allReward_R(isnan(allReward_R)) = 0;
allReward_L(isnan(allReward_L)) = 0;
allChoice_R = double(allChoices == 1);
allChoice_L = double(allChoices == -1);

allRewards = zeros(1,length(allChoices));
allRewards(logical(allReward_R)) = 1;
allRewards(logical(allReward_L)) = -1;

allITIs = [sessionData(responseInds).trialEnd] - [sessionData(responseInds).CSon];
allProbsL = [sessionData(responseInds).rewardProbL];
allProbsR = [sessionData(responseInds).rewardProbR];

%% 
figure
%suptitle(filename)
suptitle('CG15 post')

%%
rMag = 1;
nrMag = rMag/2;

% trial plot
subplot(2,6,[1:6]); hold on

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
        plot([i],0,'ko','markersize',4,'linewidth',2)
    end
%     if any(i == blockSwitch)
%         plot([i i],[-1*rMag rMag],'--','linewidth',1,'Color',[30 144 255]./255)
%     end
end

% time plot
subplot(2,6,[7:12]); hold on
xlabel('Time (min)')
for i = 1:length(sessionData)
    currTime = (sessionData(i).CSon - sessionData(1).CSon)/1000/60; %convert to min
    if strcmp(sessionData(i).trialType,'CSplus')
        if ~isnan(sessionData(i).rewardR)
            if sessionData(i).rewardR == 1 % R side rewarded
                plot([currTime currTime],[0 rMag],'k')
            else
                plot([currTime currTime],[0 nrMag],'k') % R side not rewarded
            end
        elseif ~isnan(sessionData(i).rewardL)
            if sessionData(i).rewardL == 1 % L side rewarded
                plot([currTime currTime],[-1*rMag 0],'k')
            else
                plot([currTime currTime],[-1*nrMag 0],'k')
            end
            
        else % CSplus trial but no rewardL or rewardR
            plot([currTime currTime],[-rMag rMag],'r')
        end
    else % CS minus trial
        plot([currTime],0,'ko','markersize',4,'linewidth',2)
    end
%     if any(i == blockSwitch)
%         plot([currTime currTime],[-1*rMag rMag],'--','linewidth',1,'Color',[30 144 255]./255)
%     end
end
xlim([0 currTime]);


%%
set(gcf, 'Position', get(0,'Screensize'))

savepath = [root animalName sep sessionFolder sep  'figures' sep];
if isempty(dir(savepath))
    mkdir(savepath)
end

if saveFigFlag == 1
    saveFigurePDF(gcf,[savepath sep filename '_rawBehavior'])
end

