function [choiceRwdFractions] = choiceRwdFractions_opMD(xlFile, sheet)

[root, sep] = currComputer();
[~, dayList, ~] = xlsread(xlFile, sheet);
choiceFractions = [];
rwdFractions = [];
choiceRwdFractions = [];

for i = 1: length(dayList)
    sessionName = dayList{i};
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
    
    blockSwitchL = [blockSwitchL length(behSessionData)];
    tempBlockSwitchL = blockSwitchL;
    for j = 2:length(blockSwitchL)
        subVal = sum(omitInds(tempBlockSwitchL(j-1):tempBlockSwitchL(j)));
        blockSwitchL(j:end) = blockSwitchL(j:end) - subVal;
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
    allRewards(logical(allReward_L)) = 1;
    
    
    for j = 2:length(blockSwitchL)
        choiceFractionTemp = sum(allChoice_L((blockSwitchL(j-1)+1):blockSwitchL(j))) / (blockSwitchL(j)-blockSwitchL(j-1));
        choiceFractions = [choiceFractions; choiceFractionTemp]; 
        rwdFractionTemp = sum(allReward_L((blockSwitchL(j-1)+1):blockSwitchL(j))) / sum(allRewards((blockSwitchL(j-1)+1):blockSwitchL(j)));
        rwdFractions = [rwdFractions; rwdFractionTemp];
    end
end

choiceRwdFractions = horzcat(rwdFractions, choiceFractions);