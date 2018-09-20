 function [probChoice_L, sem] = ITIswitch_opMD(xlFile, animal, category)

[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

combinedITIs = [];
combinedChoice_L = [];

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
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
    
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    allChoice_R = double(allChoices == 1);
    allChoice_L = double(allChoices == -1);

    ITIs = [behSessionData(responseInds).trialEnd] - [behSessionData(responseInds).CSon];
    
    combinedITIs = [combinedITIs ITIs(1:end-1)];
    combinedChoice_L = [combinedChoice_L allChoice_L];
    
end

logSteps = exp(linspace(log(2.5),log(30),11))*1000;
binnedITIinds = discretize(combinedITIs,logSteps);

for i = 1:11
    probChoice_L(i) = mean(combinedChoice_L(binnedITIinds == i));
    sem(i) = std(combinedChoice_L(binnedITIinds == i)) / sqrt(length(combinedChoice_L(binnedITIinds == i)));
end

figure;
errorbar(logSteps(1:10), probChoice_L(1:10), sem(1:10))

