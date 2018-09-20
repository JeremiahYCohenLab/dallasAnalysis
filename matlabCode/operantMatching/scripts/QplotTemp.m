figure; hold on; plot(CG15preS_dual.dual_tForget.Q_s(:,1), '-b')
plot(CG15preS_dual.dual_tForget.Q_l(:,1), '-c')
plot(CG15preS_dual.dual_tForget.Q_s(:,2), '-r')
plot(CG15preS_dual.dual_tForget.Q_l(:,2), 'Color', [1 0 0.5])


sessionName = ['mCG16d20170605']
[root, sep] = currComputer();

[animalName, date] = strtok(sessionName, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];

if isstrprop(sessionName(end), 'alpha')
    sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
else
    sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
end

if exist(sessionDataPath,'file')
    load(sessionDataPath)
else
    [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName);
end

%% Break session down into CS+ trials where animal responded

responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
omitInds = isnan([behSessionData.rewardTime]); 

origBlockSwitch = blockSwitch;
tempBlockSwitch = blockSwitch;
for i = 2:length(blockSwitch)
    subVal = sum(omitInds(tempBlockSwitch(i-1):tempBlockSwitch(i)));
    blockSwitch(i:end) = blockSwitch(i:end) - subVal;
end

allReward_R = [behSessionData(responseInds).rewardR]; 
allReward_L = [behSessionData(responseInds).rewardL]; 
rewProb_R = [behSessionData(responseInds).rewardProbR]; 
rewProb_L = [behSessionData(responseInds).rewardProbL]; 
allChoices = NaN(1,length(behSessionData(responseInds)));
allChoices(~isnan(allReward_R)) = 1;
allChoices(~isnan(allReward_L)) = -1;

allReward_R(isnan(allReward_R)) = 0;
allReward_L(isnan(allReward_L)) = 0;
allChoice_R = double(allChoices == 1);
allChoice_L = double(allChoices == -1);

allRewards = zeros(1,length(allChoices));
allRewards(logical(allReward_R)) = 1;
allRewards(logical(allReward_L)) = -1;

allITIs = [behSessionData(responseInds).trialEnd] - [behSessionData(responseInds).CSon];
allProbsL = [behSessionData(responseInds).rewardProbL];
allProbsR = [behSessionData(responseInds).rewardProbR];
choiceTimes = [behSessionData(responseInds).rewardTime];

figure; hold on;
for i=1:length(allChoices)
    if allChoices(i) == 1
        if allRewards(i) == 1
            plot([choiceTimes(i) choiceTimes(i)], [0.8 1], '-c')
        else
            plot([choiceTimes(i) choiceTimes(i)], [0.8 0.9], '-c')
        end
    elseif allChoices(i) == -1
        if allRewards(i) == -1
            plot([choiceTimes(i) choiceTimes(i)], [0.8 1], '-m')
        else
            plot([choiceTimes(i) choiceTimes(i)], [0.8 0.9], '-m')
        end
    end
end

for i=1:length(allChoices)
    if allChoices(i) == 1
        plot([i i], [0.8 0.9], '-b')
    elseif allChoices(i) == -1
        plot([i i], [0.8 0.9], '-r')
    end
end

if blockSwitch(end) == length(allChoices)
    blockSwitch = blockSwitch(1:end-1);
end
normKern = normpdf(-15:15,0,4);
normKern = normKern / sum(normKern);
xVals = (1:(length(normKern) + length(allChoices) - 1)) - round(length(normKern)/2);
figure; hold on; yyaxis left; plot(xVals, conv(allChoices,normKern)/max(conv(allChoices,normKern)),'k','linewidth',2);
yyaxis right; plot(dynBel_CG09_0308.basePrior.probChoice(:,1), '-b', 'linewidth', 2)
yyaxis right; hold on; plot(CG09_0308_dual.fourParams_twoLearnRates_alphaForget.probChoice(:,1), '-g', 'linewidth', 2)