function [allRewards, allChoices] = opMatchingTest(filename, saveFigFlag)

if nargin < 2
    saveFigFlag = 1;
end

% % Add current folder to working directory
% currPath = pwd;
% tmp = strfind(currPath,'\');
% addpath(genpath(currPath(1:tmp(end)))); 

[root, sep] = currComputer_operantMatching();

[animalName, date] = strtok(filename, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];
sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep filename '_sessionData.mat'];
if exist(sessionDataPath,'file')
    load(sessionDataPath)
else
    [sessionData, blockSwitch] = generateSessionData_operantMatching(filename);
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

allRewardz = ones(1,length(allChoices))*2;
allRewardz(logical(allReward_R)) = 3;
allRewardz(logical(allReward_L)) = 1;

%%

allITIs = [sessionData(responseInds).trialEnd] - [sessionData(responseInds).CSon];

% n = length(sessionData) - 1;
% allITIs = [sessionData(1:n).trialEnd] - [sessionData(1:n).CSon];
% allITIs = [NaN, allITIs];
% 
% consecResponseInds = [];
% for i = 1:(length(responseInds)- 1)
%     if responseInds(i+1) > (responseInds(i) + 1)
%         allITIs(i+1) = allITIs(i) + allITIs(i+1);
%         allITIs(i) = NaN;
%     end
% end
%consecResponseInds = consecResponseInd(~isnan


%index consecRewards and Choices

%look at time as a coregressor, so reward/time, time as a period of no
%choice or reward
%limit analysis to consecutive responses

%% 
figure
suptitle(filename)


%% linear regression model
tMax = 10;
rwdMatx = [];
for i = 1:tMax
    rwdMatx(i,:) = [NaN(1,i) allRewardz(1:end-i)];
end

choiceMatx = [];
for i = 1:tMax
    choiceMatx(i,:) = [NaN(1,i) allChoices(1:end-i)];
end

itiMatx = [];
for i = 1:tMax
    itiMatx(i,:) = [NaN(1,i) allITIs(1:end-i)];
end

itiXrwd = itiMatx.*rwdMatx;
itiXchoice = itiMatx.*choiceMatx;
itiXchoiceXrwd = itiXrwd.*choiceMatx;

subplot(2,3,1); hold on

lm = fitlm([rwdMatx]', allChoices,'linear'); rsq{1} = num2str(round(lm.Rsquared.Adjusted*100)/100);
relevInds = 2:tMax+1;
coefVals = lm.Coefficients.Estimate(relevInds);
CIbands = coefCI(lm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:tMax,coefVals,errorL,errorU,'k','linewidth',2)

lm = fitlm([rwdMatx; choiceMatx]', allChoices, 'linear'); rsq{2} = num2str(round(lm.Rsquared.Adjusted*100)/100);
relevInds = 2:tMax+1;
coefVals = lm.Coefficients.Estimate(relevInds);
CIbands = coefCI(lm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[189,183,107]/255,'linewidth',2)
xlabel('Reward n Trials Back')
ylabel('\beta Coefficient')
legend(['Rwd History Only. R^2: ' rsq{1}],['Rwd x Choice History. R^2: ' rsq{2}],'location','best')
xlim([0.5 tMax+0.5])


subplot(2,3,2); hold on

lm = fitlm([choiceMatx]', allChoices, 'linear'); rsq{3} = num2str(round(lm.Rsquared.Adjusted*100)/100);
relevInds = 2:tMax+1;
coefVals = lm.Coefficients.Estimate(relevInds);
CIbands = coefCI(lm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:tMax,coefVals,errorL,errorU,'Color',[34,139,34]/255,'linewidth',2)

lm = fitlm([rwdMatx; choiceMatx]', allChoices, 'linear');
relevInds = tMax+2:length(lm.Coefficients.Estimate);
coefVals = lm.Coefficients.Estimate(relevInds);
CIbands = coefCI(lm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[189,183,107]/255,'linewidth',2)
xlabel('Choice n Trials Back')
ylabel('\beta Coefficient')
legend(['Choice History Only. R^2: ' rsq{3}],['Rwd x Choice History. R^2: ' rsq{2}],'location','best');
xlim([0.5 tMax+0.5])


subplot(2,3,3); hold on

lm = fitlm([itiMatx]', allChoices,'linear'); rsq{1} = num2str(round(lm.Rsquared.Adjusted*100)/100);
relevInds = 2:tMax+1;
coefVals = lm.Coefficients.Estimate(relevInds);
CIbands = coefCI(lm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:tMax,coefVals,errorL,errorU,'k','linewidth',2)

lm = fitlm([rwdMatx; itiMatx]', allChoices, 'linear'); rsq{2} = num2str(round(lm.Rsquared.Adjusted*100)/100);
relevInds = 2:tMax+1;
coefVals = lm.Coefficients.Estimate(relevInds);
CIbands = coefCI(lm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[189,183,107]/255,'linewidth',2)
xlabel('Reward n Trials Back')
ylabel('\beta Coefficient')
legend(['ITI History Only. R^2: ' rsq{1}],['ITI x Reward History. R^2: ' rsq{2}],'location','best')
xlim([0.5 tMax+0.5])


subplot(2,3,4); hold on

lm = fitlm([itiXrwd]', allChoices,'linear'); rsq{1} = num2str(round(lm.Rsquared.Adjusted*100)/100);
relevInds = 2:tMax+1;
coefVals = lm.Coefficients.Estimate(relevInds);
CIbands = coefCI(lm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:tMax,coefVals,errorL,errorU,'k','linewidth',2)
xlabel('Reward n Trials Back')
ylabel('\beta Coefficient')
legend(['ITIxRwd History. R^2: ' rsq{1}],'location','best')
xlim([0.5 tMax+0.5])

subplot(2,3,5); hold on

lm = fitlm([itiXchoice]', allChoices,'linear'); rsq{1} = num2str(round(lm.Rsquared.Adjusted*100)/100);
relevInds = 2:tMax+1;
coefVals = lm.Coefficients.Estimate(relevInds);
CIbands = coefCI(lm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:tMax,coefVals,errorL,errorU,'k','linewidth',2)
xlabel('Reward n Trials Back')
ylabel('\beta Coefficient')
legend(['ITIxChoice History. R^2: ' rsq{1}],'location','best')
xlim([0.5 tMax+0.5])

subplot(2,3,6); hold on

lm = fitlm([itiXchoiceXrwd]', allChoices,'linear'); rsq{1} = num2str(round(lm.Rsquared.Adjusted*100)/100);
relevInds = 2:tMax+1;
coefVals = lm.Coefficients.Estimate(relevInds);
CIbands = coefCI(lm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:tMax,coefVals,errorL,errorU,'k','linewidth',2)
xlabel('Reward n Trials Back')
ylabel('\beta Coefficient')
legend(['ITIxRwdxChoice History. R^2: ' rsq{1}],'location','best')
xlim([0.5 tMax+0.5])





% tbl = table(allRewardz', allChoices', allITIs', 'VariableNames', {'allRewards', 'allChoices', 'allITIs'});
% lm = fitlm(tbl, 'interactions', 'ResponseVar', 'allChoices', 'PredictorVars', {'allRewards', 'allITIs'});
% plotInteraction(lm,allRewards,allChoices);


%make sure it's comparing the previous ITI
%find a way to plot it



