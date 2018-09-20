function kilosort_opMD(sessionName, recType, filename)

% Path
[root,sep] = currComputer();

[animalName] = strtok(sessionName, 'd');
animalName = animalName(2:end);

if isstrprop(sessionName(end), 'alpha')
    fPath = [root animalName sep sessionName(1:end-1) sep 'ephys' sep 'session ' sessionName(end) sep];
    savePath = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
else
    fPath = [root animalName sep sessionName sep 'ephys' sep recType sep];
    savePath = [root animalName sep sessionName sep 'sorted' sep recType sep];
end

addpath(genpath('C:\Users\cooper_PC\Desktop\GitHub repositories\cooperAnalysis\kilosort')) % path to kilosort folder
addpath(genpath('C:\Users\cooper_PC\Desktop\GitHub repositories\cooperAnalysis\python\npy-matlab-master')) % path to npy-matlab scripts

pathToYourConfigFile = 'C:\Users\cooper_PC\Desktop\GitHub repositories\cooperAnalysis\matlabCode\operantMatching\kilosort_opMD'; % take from Github folder and put it somewhere else (together with the master_file)
run(fullfile(pathToYourConfigFile, 'StandardConfig_opMD.m'))

tic; % start timer
%
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
end
%
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this)
%     rez = merge_posthoc2(rez);

% save matlab results file
%use make path statement here
save(fullfile(savePath, 'rez.mat'), 'rez', '-v7.3');

% save python results file for Phy
rezToPhy(rez, ops.root);

% remove temporary file
delete(ops.fproc);
%%
