function [combinedStruct, mdl] = combineStructQ_opM(xlFile, animal, category)

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
    
    behDataPath = [root animalName sep sessionFolder sep 'behavior' sep sessionName '.asc']; 

    if exist(sessionDataPath,'file')
        load(sessionDataPath)
    else
        [sessionData, blockSwitch, blockProbs] = generateSessionData_behav_operantMatching(behDataPath);
    end
    
    combinedStruct = [combinedStruct sessionData];
    combinedBlockSwitch = [combinedBlockSwitch blockSwitch];
    
end

sessionData = combinedStruct;
blockSwitch = combinedBlockSwitch;

if strfind(category, 'preS')
    date = 'd20182222';
elseif strfind(category, 'postS')
    date = 'd20183333';
elseif strfind(category, 'pre')
    date = 'd20184444';
elseif strfind(category, 'post')
    date = 'd20185555';
end


savePath = [root animalName sep 'm' animalName date sep 'sorted' sep 'session' sep];
if isempty(dir(savePath))
    mkdir(savePath)
end

save([savePath 'm' animalName date '_sessionData_behav.mat'], 'sessionData', 'blockSwitch');




