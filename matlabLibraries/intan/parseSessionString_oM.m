function pathData = parseSessionString_oM(fileOrFolder, root, sep)
% parseSessionString_oM    Parses input string to generate corresponding pathData outputs
%   INPUTS
%       fileOrFolder: sessionName or name of .asc file
%           e.g.: 'mBB041d20161006' or 'mBB041d20161006.asc'
%       root: root folder
%           e.g.: 'G:\'
%       sep: separator
%           e.g.: '\' or '/'
%   OUTPUTS
%       pathData
%           Structure with sessionFolder, sortedFolder, etc...

if contains(fileOrFolder,'.asc') % input is .asc file
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
    fileInd = contains({allFiles.name},'.asc');
    behavioralDataPath = [filepath allFiles(fileInd).name];
    suptitleName = allFiles(fileInd).name(1:end-4);
    saveFigName = sessionFolder(~(sessionFolder==sep));
end

if isstrprop(suptitleName(end), 'alpha')
    sortedFolderLocation = [root animalName sep sessionFolder sep 'sorted' sep 'session ' suptitleName(end) sep];
    pathData.saveFigFolder = [root animalName sep sessionFolder sep 'figures' sep 'session ' suptitleName(end) sep];
else
    sortedFolderLocation = [root animalName sep sessionFolder sep 'sorted' sep];
    pathData.saveFigFolder = [root animalName sep sessionFolder sep 'figures' sep];
end

% append path information
pathData.suptitleName = suptitleName;
pathData.sessionFolder = sessionFolder;
pathData.sortedFolder = sortedFolderLocation;
pathData.animalName = animalName;
pathData.saveFigName = saveFigName;
pathData.baseFolder = [root animalName sep sessionFolder sep];
pathData.behavioralDataPath = behavioralDataPath;
pathData.ephysPath = [root animalName sep sessionFolder sep 'ephys' sep];