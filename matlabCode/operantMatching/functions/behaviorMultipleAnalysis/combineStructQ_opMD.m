function [combinedStruct, mdl] = combineStructQ_opMD(xlFile, animal, category)

[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
combinedStruct = []; 
combinedBlockSwitch = []; 

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
    
    combinedStruct = [combinedStruct behSessionData];
    combinedBlockSwitch = [combinedBlockSwitch blockSwitch];
    
end

behSessionData = combinedStruct;
blockSwitch = combinedBlockSwitch;

if strfind(category, 'preS')
    date = 'd20177777';
elseif strfind(category, 'pre')
    date = 'd20178888';
elseif strfind(category, 'post')
    date = 'd20179999';
elseif strfind(category, 'effect')
    date = 'd20176666';
end


savePath = [root animalName sep 'm' animalName date sep 'sorted' sep 'session' sep];
if isempty(dir(savePath))
    mkdir(savePath)
end

save([savePath 'm' animalName date '_sessionData_behav.mat'], 'behSessionData', 'blockSwitch');

mdl = qLearning_fit2LR(['m' animalName date '.asc']);
alphaNPE = mdl.threeParams_twoLearnRates.bestParams(1);
alphaPPE = mdl.threeParams_twoLearnRates.bestParams(2);
beta = mdl.threeParams_twoLearnRates.bestParams(3);




