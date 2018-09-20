function trialSpikeRasters_opMD(sessionName, saveFigFlag)

if nargin < 2
    saveFigFlag = 1
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

sortedFolder = dir(sortedFolderLocation);
    
if any(~cellfun(@isempty,strfind({sortedFolder.name},'_behav.mat'))) == 1
    sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_behav.mat'));
    load([sortedFolderLocation sortedFolder(sessionDataInd).name],'blockSwitch', 'blockSwitchL', 'blockSwitchR', 'behSessionData');
else  
    [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName);
end

if any(~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'))) == 1
    sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'));
    load([sortedFolderLocation sortedFolder(sessionDataInd).name])
else
    [sessionData] = generateSessionData_nL_operantMatching(sessionName);
end

if isempty(dir(savepath))
    mkdir(savepath)
end

tb = 1.5;
tf = 8;
time = -1000*tb:1000*tf;

omitInds = isnan([sessionData.rewardTime]);
tempBlockSwitch = blockSwitch;
blockSwitch_CSminCorrected = blockSwitch;
for i = 2:length(blockSwitch_CSminCorrected)
    subVal = sum(omitInds(tempBlockSwitch(i-1):tempBlockSwitch(i)));
    blockSwitch_CSminCorrected(i:end) = blockSwitch_CSminCorrected(i:end) - subVal;
end



%% Sort all spikes into a raster-able matrix
spikeFields = fields(sessionData);
clust = find(~cellfun(@isempty,strfind(spikeFields,'SS')) & ~cellfun(@isempty,strfind(spikeFields,'TT')));

allTrial_spike = {};
for k = 1:length(sessionData)
    for i = 1:length(clust)
        if k == 1
            prevTrial_spike = [];
            currTrial_lick = [];
        else
            prevTrial_spikeInd = sessionData(k-1).(spikeFields{clust(i)}) > sessionData(k-1).trialEnd-tb*1000;
            prevTrial_spike = sessionData(k-1).(spikeFields{clust(i)})(prevTrial_spikeInd) - sessionData(k).CSon;
        end
        
        currTrial_spikeInd = sessionData(k).(spikeFields{clust(i)}) < sessionData(k).CSon+tf*1000;
        currTrial_spike = sessionData(k).(spikeFields{clust(i)})(currTrial_spikeInd) - sessionData(k).CSon;
        
        allTrial_spike{i,k} = [prevTrial_spike currTrial_spike];

    end
    if ~isnan(sessionData(k).licksL)
        currTrial_lickInd = sessionData(k).licksL < sessionData(k).CSon + tf*1000;
        currTrial_lick = sessionData(k).licksL(currTrial_lickInd) - sessionData(k).CSon;
    elseif ~isnan(sessionData(k).licksR)
        currTrial_lickInd = sessionData(k).licksR < sessionData(k).CSon + tf*1000;
        currTrial_lick = sessionData(k).licksR(currTrial_lickInd) - sessionData(k).CSon;  
    else
        currTrial_lick = 0;
    end
    allTrial_lick{k} = [currTrial_lick];
end

% sometimes no licks/spikes are considered 1x0 and sometimes they are []
% plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
allTrial_spike(cellfun(@isempty,allTrial_spike)) = {zeros(1,0)}; 

csXinds = [0 0 500 500];
lickStruct.Color = [0 1 1];
lickStruct.LineWidth = 3;
spikeStruct.LineWidth = 1.5;

for i = 1:length(clust)
    rawSpikes = figure;
    titleTmp = strrep(spikeFields{clust(i)},'_','\_');
    suptitle([sessionName ': ' titleTmp]);
    plot(repmat([-5000 10000],length(blockSwitch),1)', [blockSwitch; blockSwitch],'b'); hold on
    plotSpikeRaster(allTrial_lick(1,:),'PlotType','vertline','LineFormat', lickStruct); hold on
    plotSpikeRaster(allTrial_spike(i,:),'PlotType','vertline','LineFormat', spikeStruct); hold on;
    csYinds = [0.5 1.5 1.5 0.5];
    for j = 1:length(allTrial_spike)
        if regexp(behSessionData(j).trialType,'CSplus')
            csPlus(j) = fill(csXinds,csYinds,'b'); set(csPlus(j), 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
            if behSessionData(j).rewardL == 1 || sessionData(j).rewardR == 1
                rwdTmp = (behSessionData(j).rewardTime - behSessionData(j).CSon);
                rwdTimeInds = [rwdTmp rwdTmp (rwdTmp+50) (rwdTmp+50)];
                rwdTime(j) = fill(rwdTimeInds,csYinds,[0.5 0 0.8]); set(rwdTime(j), 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
            end
        else
            csMinus(j) = fill(csXinds,csYinds,'r'); set(csMinus(j), 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
        end
        csYinds = csYinds + 1; 
    end
    h = zeros(5, 1);
    h(1) = plot(NaN,NaN,'sk', 'MarkerFaceColor', 'k');
    h(2) = plot(NaN,NaN,'s', 'Color', [0.5 0.5 1], 'MarkerFaceColor',[0.5 0.5 1]);
    h(3) = plot(NaN,NaN,'s', 'Color', [1 0.5 0.5], 'MarkerFaceColor',[1 0.5 0.5]);
    h(4) = plot(NaN,NaN,'s', 'Color', [0 1 1], 'MarkerFaceColor',[0 1 1]);
    h(5) = plot(NaN,NaN,'s', 'Color', [0.5 0 0.8], 'MarkerFaceColor',[0.5 0 0.8]);
    legend(h, 'spikes','CS plus','CS minus', 'licks', 'reward');    
    if saveFigFlag == 1 
        saveFigurePDF(rawSpikes,[savepath sep sessionName '_' spikeFields{clust(i)}, '_raw'])
    end
end

end