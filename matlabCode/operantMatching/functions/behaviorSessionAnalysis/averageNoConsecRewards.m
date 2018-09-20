function [allRewards, allChoices] = averageNoConsecRewards(filename, saveFigFlag)


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
changeChoice = [false abs(diff(allChoices)) > 0];
changeHistogram = [];
changeHistogram_LtoR = [];
changeHistogram_RtoL = [];
for i = find(changeChoice == 1)
    if allChoices(i) == 1 % if a right lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allChoices(i-goBack) == -1 && allRewards(i-goBack) == 0 % if previous trial was a L lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram = [changeHistogram temp];
        changeHistogram_LtoR = [changeHistogram_LtoR temp];
    elseif allChoices(i) == -1 %if a left lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allChoices(i-goBack) == 1 && allRewards(i-goBack) == 0 % if previous trial was a R lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram = [changeHistogram temp];
        changeHistogram_RtoL = [changeHistogram_RtoL temp];
    end
end


figure
