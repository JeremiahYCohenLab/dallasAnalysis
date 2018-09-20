function [allRewards, allChoices] = opMatchingPerseveration(filename)

% % Add current folder to working directory
% currPath = pwd;
% tmp = strfind(currPath,'\');
% addpath(genpath(currPath(1:tmp(end)))); 

[root, sep] = currComputer();

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


%%

normKern = normpdf(-4:4,0,4);
normKern = normKern / sum(normKern);
halfKern = normKern(round(length(normKern)/2):end);
xVals = (1:(length(halfKern) + length(allChoices) - 1)) - round(length(halfKern)/2);
meowChoices = conv(allChoices,halfKern)/max(conv(allChoices,halfKern));
meowRewards = conv(allRewards,halfKern)/max(conv(allRewards,halfKern));

meowRewChoice = vertcat(meowRewards, meowChoices);


% k = 5;
% l = k-1;
% divLength = floor(length(allChoices)/k);
% 
% for i = 1:divLength
%     j=i*k;
%     fractionChoice_L(i) = sum(allChoice_L((j-l):j)) / (sum(allChoice_R((j-l):j)) + sum(allChoice_L((j-l):j)));
%     fractionReward_L(i) = sum(allReward_L((j-l):j)) / (sum(allReward_R((j-l):j)) + sum(allReward_L((j-l):j)));
% end
% 
% fractionReward_L(isnan(fractionReward_L)) = 0;
% fractionsChoiceReward = vertcat(fractionReward_L, fractionChoice_L);

% normKern = ones(1, 10);
% temp = conv(allChoice_L,normKern)/max(conv(allChoice_L,normKern));
% smoothedChoice_L = temp(1:length(allChoices));
% temp = conv(allChoice_R,normKern)/max(conv(allChoice_R,normKern));
% smoothedChoice_R = temp(1:length(allChoices));
% 
% temp = conv(allReward_L,normKern)/max(conv(allReward_L,normKern));
% smoothedReward_L = temp(1:length(allRewards));
% temp = conv(allReward_R,normKern)/max(conv(allReward_R,normKern));
% smoothedReward_R = temp(1:length(allRewards));
% 
% fractionChoice = smoothedChoice_L ./ (smoothedChoice_L + smoothedChoice_R);
% fractionReward = smoothedReward_L ./ (smoothedReward_L + smoothedReward_R);
% fractionReward(isnan(fractionReward)) = 0;


%%
choiceBlocks = [];
choiceBlocks_L = [];
choiceBlocks_R = [];
rewardBlocks = [];
rewardBlocks_L = [];
rewardBlocks_R = [];
blockSwitch = blockSwitch + 1;
blockSwitch(1) = 1;
lastTrial = length(allChoices);
blockSwitchAndEnd = horzcat(blockSwitch, lastTrial);
totalSwitches = length(blockSwitchAndEnd);

j = 1;
k = 1;
for i = 2:(totalSwitches - 2)
    if rewProb_L(blockSwitchAndEnd(i)) > rewProb_R(blockSwitchAndEnd(i))
        choiceBlocks(i-1, :) = allChoice_L(blockSwitch(i):(blockSwitch(i)+69));
        rewardBlocks(i-1, :) = allReward_L(blockSwitch(i):(blockSwitch(i)+69));
        choiceBlocks_L(j, :) = allChoice_L(blockSwitch(i):(blockSwitch(i)+69));
        rewardBlocks_L(j, :) = allReward_L(blockSwitch(i):(blockSwitch(i)+69));
        j = j + 1;
    else
        choiceBlocks(i-1, :) = allChoice_R(blockSwitch(i):(blockSwitch(i)+69));
        rewardBlocks(i-1, :) = allReward_R(blockSwitch(i):(blockSwitch(i)+69));
        choiceBlocks_R(k, :) = allChoice_R(blockSwitch(i):(blockSwitch(i)+69));
        rewardBlocks_R(k, :) = allReward_R(blockSwitch(i):(blockSwitch(i)+69));
        k = k + 1;
    end
end

if (blockSwitchAndEnd(totalSwitches) - blockSwitchAndEnd(totalSwitches - 1) > 68)
    if rewProb_L(blockSwitchAndEnd(totalSwitches - 1)) > rewProb_R(blockSwitchAndEnd(totalSwitches - 1))
        choiceBlocks(totalSwitches-2, :) = allChoice_L(blockSwitch(end):(blockSwitch(end) + 69));
        rewardBlocks(totalSwitches-2, :) = allReward_L(blockSwitch(end):(blockSwitch(end) + 69));
        choiceBlocks_L(k, :) = allChoice_L(blockSwitch(end):(blockSwitch(end) + 69));
        rewardBlocks_L(k, :) = allReward_L(blockSwitch(i):(blockSwitch(i)+69));
    else
        choiceBlocks(totalSwitches-2, :) = allChoice_R(blockSwitch(end):(blockSwitch(end) + 69));
        rewardBlocks(totalSwitches-2, :) = allReward_R(blockSwitch(end):(blockSwitch(end) + 69));
        choiceBlocks_R(k, :) = allChoice_R(blockSwitch(end):(blockSwitch(end) + 69));
        rewardBlocks_R(k, :) = allReward_R(blockSwitch(i):(blockSwitch(i)+69));
    end
end


tempz = choiceBlocks - rewardBlocks;
height = size(tempz, 1);
boxKern = ones(1, 5);

% height2 = size(choiceBlocks_L,1);
% 
% for i = 1:height2
%     temp = conv(choiceBlocks_L(i,:),normKern)/max(conv(choiceBlocks_L(i,:),normKern));
%     meowChoices(i,:) = temp(1:70);
%     temp = conv(rewardBlocks_L(i,:),normKern)/max(conv(rewardBlocks_L(i,:),normKern));
%     meowRewards(i,:) = temp(1:70);
% end
% 
% meowDiff = meowChoices - meowRewards;
% meowSEM = std(meowDiff, 1) / sqrt(length(meowDiff));
% meowAverage = mean(meowDiff, 1);
% aaa_meow = vertcat(meowAverage, meowSEM);

for i = 1:height
    temp = conv(tempz(i,:),boxKern)/max(conv(tempz(i,:),boxKern));
    diffChoiceRewards(i,:) = temp(1:70);
end

sem = std(diffChoiceRewards, 1) / sqrt(length(diffChoiceRewards));
diffChoiceRewardsAve = mean(diffChoiceRewards, 1);
diffChoiceRewardsandSEM = vertcat(diffChoiceRewardsAve, sem);

tempz_L = choiceBlocks_L - rewardBlocks_L;
height = size(tempz_L, 1);

for i = 1:height
    temp = conv(tempz_L(i,:),boxKern)/max(conv(tempz_L(i,:),boxKern));
    diffChoiceRewards_L(i,:) = temp(1:70);
end

sem_L = std(diffChoiceRewards_L, 1) / sqrt(length(diffChoiceRewards_L));
diffChoiceRewardsAve_L = mean(diffChoiceRewards_L, 1);
diffChoiceRewardsandSEM_L = vertcat(diffChoiceRewardsAve_L, sem);

tempz_R = choiceBlocks_R - rewardBlocks_R;
height = size(tempz_R, 1);

for i = 1:height
    temp = conv(tempz_R(i,:),boxKern)/max(conv(tempz_R(i,:),boxKern));
    diffChoiceRewards_R(i,:) = temp(1:70);
end

sem_R = std(diffChoiceRewards_R, 1) / sqrt(length(diffChoiceRewards_R));
diffChoiceRewardsAve_R = mean(diffChoiceRewards_R, 1);
diffChoiceRewardsandSEM_R = vertcat(diffChoiceRewardsAve_R, sem_R);



figure