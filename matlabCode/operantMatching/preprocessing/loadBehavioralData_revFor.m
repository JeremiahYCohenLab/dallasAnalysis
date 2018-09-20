function [behSessionData, blockSwitch, outputPathStruct] = loadBehavioralData_revFor(fileOrFolder)

[root, sep] = currComputer();

if ~isempty(strfind(fileOrFolder,'.asc')) % input is .asc file
    filename = fileOrFolder;
    [animalName, date] = strtok(filename, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];
    behavioralDataPath = [root animalName sep sessionFolder sep 'behavior' sep filename];
    suptitleName = filename(1:strfind(filename,'.asc')-1);
    saveFigName = suptitleName;
else % input is the folder
    sessionFolder = fileOrFolder;
    animalName = strtok(sessionFolder, 'd');
    animalName = animalName(2:end);
    filepath = [root animalName sep sessionFolder sep 'behavior' sep];
    allFiles = dir(filepath);
    fileInd = ~cellfun(@isempty,strfind({allFiles.name},'.asc'));
    behavioralDataPath = [filepath allFiles(fileInd).name];
    suptitleName = allFiles(fileInd).name;
    saveFigName = sessionFolder(~(sessionFolder==sep));
end

if isstrprop(suptitleName(end), 'alpha')
    sortedFolderLocation = [root animalName sep sessionFolder sep 'sorted' sep 'session ' suptitleName(end) sep];
else
    sortedFolderLocation = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep];
end
sortedFolder = dir(sortedFolderLocation);
sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_behav.mat')) & ~cellfun(@isempty,strfind({sortedFolder.name},suptitleName)); 
if any(sessionDataInd) % check if there is a file with suptitleName prefix and _behav.mat suffix
    load([sortedFolderLocation sortedFolder(sessionDataInd).name])
    behSessionData = sessionData;
else
    [behSessionData, blockSwitch, ~] = generateSessionData_behav_operantMatching(behavioralDataPath);
end

outputPathStruct.suptitleName = suptitleName;
outputPathStruct.sessionFolder = sessionFolder;
outputPathStruct.animalName = animalName;
outputPathStruct.saveFigName = saveFigName;