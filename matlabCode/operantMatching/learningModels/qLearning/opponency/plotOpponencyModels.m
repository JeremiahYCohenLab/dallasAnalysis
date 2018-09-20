function plotOpponencyModels(sessionName)


[mdl] = qLearning_fit2LR(strcat(sessionName,'.asc'));

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

if blockSwitch(end) == length(allChoices)
    blockSwitch = blockSwitch(1:end-1);
end


%% 
figure
suptitle(sessionName)


normKern = normpdf(-15:15,0,4);
normKern = normKern / sum(normKern);
xVals = (1:(length(normKern) + length(allChoices) - 1)) - round(length(normKern)/2);

subplot(4,6,[1:6]); hold on
yyaxis left; plot(xVals, conv(allChoices,normKern)/max(conv(allChoices,normKern)),'k','linewidth',2);
ylabel('<-- Left       Right -->')
yyaxis right; plot(mdl.twoParams.probChoice(:,1), '-', 'Color', [0 1 1], 'linewidth', 2)
ylabel('prob R choice')
legend('Choices','2 Param Q')
xlim([1 length(allChoice_R)])


subplot(4,6,[7:12]); hold on
set(gca,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
yyaxis left; plot(xVals, conv(allChoices,normKern)/max(conv(allChoices,normKern)),'k','linewidth',2);
yyaxis right; plot(mdl.fourParams_twoLearnRates_alphaForget.probChoice(:,1), '-g', 'linewidth', 2)
legend('Choices','4 param Q')
xlim([1 length(allChoice_R)])



subplot(4,6,[13:18]); hold on
yyaxis left; plot(xVals, conv(allChoices,normKern)/max(conv(allChoices,normKern)),'k','linewidth',2);
yyaxis right; plot(mdl.fourParams_opponency.probChoice(:,1), '-b', 'linewidth', 2)
legend('Choices','4 Param PAR Q')
xlim([1 length(allChoice_R)])


subplot(4,6,[19:24]); hold on
yyaxis left; plot(xVals, conv(allChoices,normKern)/max(conv(allChoices,normKern)),'k','linewidth',2);
yyaxis right; plot(mdl.fiveParams_opponency.probChoice(:,1), '-','Color', [0.7 0 1], 'linewidth', 2)
xlabel('Trials')
legend('Choices','5 Param PAR Q')
xlim([1 length(allChoice_R)])


