function [spikeStruct] = spikeProps_opMD(sessionName, cellName, saveFigFlag)

if nargin < 3
    saveFigFlag = 1;
end

% Path
[root,sep] = currComputer();

[animalName] = strtok(sessionName, 'd');
animalName = animalName(2:end);

if isstrprop(sessionName(end), 'alpha')
    sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
    savepath = [root animalName sep sessionName(1:end-1) sep  'figures' sep 'session ' sessionName(end) sep];
else
    sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
    savepath = [root animalName sep sessionName sep  'figures' sep];
end
if isempty(dir(savepath))
    mkdir(savepath)
end

sortedFolder = dir(sortedFolderLocation);

if any(~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'))) == 1
    sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'));
    load([sortedFolderLocation sortedFolder(sessionDataInd).name])
else
    [sessionData] = generateSessionData_nL_operantMatching(sessionName);
end


%% Sort all spikes into a raster-able matrix

%set time window for spike analyses
tb = 1.5;
tf = 8;
time = -1000*tb:1000*tf;

spikeFields = fields(sessionData);
cellInd = find(~cellfun(@isempty,strfind(spikeFields,cellName)));


allTrial_spike = {};
for k = 1:length(sessionData)
    if k == 1
        prevTrial_spike = [];
    else
        prevTrial_spikeInd = sessionData(k-1).(spikeFields{cellInd}) > sessionData(k-1).trialEnd-tb*1000;
        prevTrial_spike = sessionData(k-1).(spikeFields{cellInd})(prevTrial_spikeInd) - sessionData(k).CSon;
    end

    currTrial_spikeInd = sessionData(k).(spikeFields{cellInd}) < sessionData(k).CSon+tf*1000;
    currTrial_spike = sessionData(k).(spikeFields{cellInd})(currTrial_spikeInd) - sessionData(k).CSon;

    allTrial_spike{k} = [prevTrial_spike currTrial_spike];
end

% sometimes no spikes are considered 1x0 and sometimes they are []
% plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
allTrial_spike(cellfun(@isempty,allTrial_spike)) = {zeros(1,0)}; 


%%

smoothWin = 250;
trialBeg = tb*1000;
CSoff = tb*1000 + 500;
for j = 1:length(sessionData)
    trialDurDiff(j) = (sessionData(j).trialEnd - sessionData(j).CSon)- tf*1000;
end



allTrial_spikeMatx = NaN(length(sessionData),length(time));
for j = 1:length(allTrial_spike)
    tempSpike = allTrial_spike{j};
    tempSpike = tempSpike + tb*1000; % add this to pad time for SDF
    allTrial_spikeMatx(j,tempSpike) = 1;
    if trialDurDiff(j) < 0
        allTrial_spikeMatx(j, isnan(allTrial_spikeMatx(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
    else
        allTrial_spikeMatx(j, isnan(allTrial_spikeMatx(j,:))) = 0;
    end
    if sum(allTrial_spikeMatx(j,:)) == 0     %if there is no spike data for this trial, don't count it
        allTrial_spikeMatx(j,:) = NaN;
    end
end


for j = 1:length(allTrial_spike)
    if ~isempty(allTrial_spikeMatx(j))
        preCSspikeCount(j) = sum(allTrial_spikeMatx(j, 1:trialBeg));              %find total spikes before CS on
        postCSspikeCount(j) = sum(allTrial_spikeMatx(j, trialBeg:CSoff));    

        spikeTemp = fastsmooth(allTrial_spikeMatx(j,:)*1000, smoothWin, 3);         %smooth raw spikes to find features of spike rate
        maxFRcs(j) = max(spikeTemp(trialBeg:CSoff));
        minFRcs(j) = min(spikeTemp(trialBeg:CSoff));
        if ~isnan(maxFRcs(j))
            maxFRtime(j) = find(spikeTemp == max(spikeTemp(trialBeg:CSoff)), 1);
        else
            maxFRtime(j) = NaN;
        end

    else
        preCSspikeCount(j) = NaN;
        postCSspikeCount(j) = NaN;
        maxFRcs(j) = NaN;
        minFRcs(j) = NaN;
        maxFRtime(j) = NaN;
    end
end

responseInds = find(~isnan([sessionData.rewardTime]));
spikeStruct = struct;
spikeStruct.preCScount = preCSspikeCount(responseInds);
spikeStruct.postCScount = postCSspikeCount(responseInds);
spikeStruct.maxCSrate = maxFRcs(responseInds);
spikeStruct.minCSrate = minFRcs(responseInds);



    