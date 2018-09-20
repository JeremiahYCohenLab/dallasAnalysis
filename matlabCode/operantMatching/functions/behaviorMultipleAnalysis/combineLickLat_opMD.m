function [lickLat_L, lickLat_R] = combineLickLat_opMD(xlFile, animal, category)

[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
lickLat_L = [];
lickLat_R = [];

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
    
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allChoice_R = double(allChoices == 1);
    allChoice_L = double(allChoices == -1);
    
    lickLatTemp_L = nan(1, length(allChoices));
    lickLatTemp_R = nan(1, length(allChoices));
    
    for j = 1:length(allChoices)
        if allChoices(j) == -1
            lickLatTemp_L(j) = behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j)).CSon;
        else
            lickLatTemp_R(j) = behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j)).CSon;
        end
    end
    
    lickLat_L = [lickLat_L lickLatTemp_L(~isnan(lickLatTemp_L))];
    lickLat_R = [lickLat_R lickLatTemp_R(~isnan(lickLatTemp_R))];
end
