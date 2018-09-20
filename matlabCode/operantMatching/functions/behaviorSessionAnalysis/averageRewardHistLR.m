function [allRewards, allChoices] = averageRewardHistLR(filename, saveFigFlag)

if nargin < 2
    saveFigFlag = 1;
end

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

normKern = normpdf(-4:4,0,4);
normKern = normKern / sum(normKern);
halfKern = normKern(round(length(normKern)/2):end);
xVals_R = (1:(length(halfKern) + length(allReward_R) - 1)) - round(length(halfKern)/2);
xVals_L = (1:(length(halfKern) + length(allReward_L) - 1)) - round(length(halfKern)/2);

smoothedReward_R = conv(allReward_R,halfKern)/max(conv(allReward_R,halfKern));
smoothedReward_L = conv(allReward_L,halfKern)/max(conv(allReward_L,halfKern));

smoothedReward_R = smoothedReward_R(2:(length(smoothedReward_R) - 4));
smoothedReward_L = smoothedReward_L(2:(length(smoothedReward_L) - 4));

choices = allChoices(2:length(allChoices));

copypasta_R = vertcat(choices, smoothedReward_R);
copypasta_L = vertcat(choices, smoothedReward_L);


a = [1];
    

