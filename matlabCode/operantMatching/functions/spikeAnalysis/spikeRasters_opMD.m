function spikeRasters_opMD(sessionName, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('cellName', ['all']);
p.addParameter('saveFigFlag', 1)
p.addParameter('intanFlag',0)
p.addParameter('revForFlag',0)
p.addParameter('modelsFlag',0)
p.addParameter('timeMax', 121000)
p.addParameter('timeBins', 12)
p.addParameter('tb', 1.5);
p.addParameter('tf', 3);
p.parse(varargin{:});

cellName = p.Results.cellName;

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
    

if p.Results.intanFlag == 1
    if exist([sortedFolderLocation sessionName '_sessionData_intan.mat'],'file')
        load([sortedFolderLocation sessionName '_sessionData_intan.mat'])
    else
        [sessionData] = generateSessionData_intan_operantMatching(sessionName);
    end
else
    if exist([sortedFolderLocation sessionName '_sessionData_nL.mat'],'file')
        load([sortedFolderLocation sessionName '_sessionData_nL.mat'])
    else
        [sessionData] = generateSessionData_nL_operantMatching(sessionName);
    end
end

[s] = behAnalysisNoPlot_opMD(sessionName, 'revForFlag', p.Results.revForFlag);

if isempty(dir(savepath))
    mkdir(savepath)
end


time = -1000*p.Results.tb:1000*p.Results.tf;

omitInds = isnan([sessionData.rewardTime]);
tempBlockSwitch = s.blockSwitch;
blockSwitch_CSminCorrected = s.blockSwitch;
for i = 2:length(blockSwitch_CSminCorrected)
    subVal = sum(omitInds(tempBlockSwitch(i-1):tempBlockSwitch(i)));
    blockSwitch_CSminCorrected(i:end) = blockSwitch_CSminCorrected(i:end) - subVal;
end


%% Sort all spikes into a raster-able matrix
spikeFields = fields(sessionData);
if iscell(cellName)
    for i = 1:length(cellName)
        clust(i) = find(~cellfun(@isempty,strfind(spikeFields,cellName{i})));
    end
elseif regexp(cellName, 'all')
    if p.Results.intanFlag
        clust = find(~cellfun(@isempty,strfind(spikeFields,'C_')));
    else
        clust = find(~cellfun(@isempty,strfind(spikeFields,'SS')) & ~cellfun(@isempty,strfind(spikeFields,'TT')));
    end
else
    clust = find(~cellfun(@isempty,strfind(spikeFields,cellName)));
end
    
allTrial_spike = {};
for k = 1:length(sessionData)
    for i = 1:length(clust)
        if k == 1
            prevTrial_spike = [];
            currTrial_lick = [];
        else
            prevTrial_spikeInd = sessionData(k-1).(spikeFields{clust(i)}) > sessionData(k-1).trialEnd-p.Results.tb*1000;
            prevTrial_spike = sessionData(k-1).(spikeFields{clust(i)})(prevTrial_spikeInd) - sessionData(k).CSon;
        end
        
        currTrial_spikeInd = sessionData(k).(spikeFields{clust(i)}) < sessionData(k).CSon+p.Results.tf*1000;
        currTrial_spike = sessionData(k).(spikeFields{clust(i)})(currTrial_spikeInd) - sessionData(k).CSon;
        
        allTrial_spike{i,k} = [prevTrial_spike currTrial_spike];

    end
    if ~isnan(sessionData(k).licksL)
        currTrial_lickInd = sessionData(k).licksL < sessionData(k).CSon + p.Results.tf*1000;
        currTrial_lick = sessionData(k).licksL(currTrial_lickInd) - sessionData(k).CSon;
    elseif ~isnan(sessionData(k).licksR)
        currTrial_lickInd = sessionData(k).licksR < sessionData(k).CSon + p.Results.tf*1000;
        currTrial_lick = sessionData(k).licksR(currTrial_lickInd) - sessionData(k).CSon;  
    else
        currTrial_lick = 0;
    end
    allTrial_lick{k} = [currTrial_lick];
end

% sometimes no licks/spikes are considered 1x0 and sometimes they are []
% plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
allTrial_spike(cellfun(@isempty,allTrial_spike)) = {zeros(1,0)}; 



%% Generate indices for block switches, trial type, choice, outcomes
blockS = [];
blockL = [];
blockR = [];
s.blockSwitch = s.blockSwitch + 1;
s.blockSwitch(1) = 1; %% this is to avoid an error on the first iteration of the loop because it's trying to index at 0

if p.Results.revForFlag
    for i = 1:length(s.blockSwitch)
        [rewardProbL, rewardProbR] = strtok(s.blockProbs(i), '/');
        rewardProbL = str2double(rewardProbL); rewardProbR = str2double(rewardProbR{1}(2:end));
        if rewardProbL > rewardProbR
            if i ~= length(s.blockSwitch)
                blockL = [blockL s.blockSwitch(i):s.blockSwitch(i+1)-1];
            else
                blockL = [blockL s.blockSwitch(i):length(sessionData)];
            end
        else
            if i ~= length(s.blockSwitch)
                blockR = [blockR s.blockSwitch(i):s.blockSwitch(i+1)-1];
            else
                blockR = [blockR s.blockSwitch(i):length(sessionData)];
            end
        end
    end
else
    for i = 1:length(s.blockSwitch)
        if s.behSessionData(s.blockSwitch(i)).rewardProbL == s.behSessionData(s.blockSwitch(i)).rewardProbR   
            if i ~= length(s.blockSwitch)
                blockS = [blockS s.blockSwitch(i):s.blockSwitch(i+1)-1];
            else
                blockS = [blockS s.blockSwitch(i):length(sessionData)];
            end
        elseif s.behSessionData(s.blockSwitch(i)).rewardProbL > s.behSessionData(s.blockSwitch(i)).rewardProbR
            if i ~= length(s.blockSwitch)
                blockL = [blockL s.blockSwitch(i):s.blockSwitch(i+1)-1];
            else
                blockL = [blockL s.blockSwitch(i):length(sessionData)];
            end
        else
            if i ~= length(s.blockSwitch)
                blockR = [blockR s.blockSwitch(i):s.blockSwitch(i+1)-1];
            else
                blockR = [blockR s.blockSwitch(i):length(sessionData)];
            end
        end
    end
end


% CSplus trials
CSplus_Inds = strcmp({sessionData.trialType},'CSplus');

% CSminus trials
CSminus_Inds = strcmp({sessionData.trialType},'CSminus');

% Lick R trials
Rlicks = s.lickSide == 1;
[~, lickR_Inds] = sort(s.lickLat);
lickR_Inds = lickR_Inds(Rlicks(lickR_Inds));

% Lick L trials
Llicks = s.lickSide == -1;
[~, lickL_Inds] = sort(s.lickLat);
lickL_Inds = lickL_Inds(Llicks(lickL_Inds));

% Rewarded trials
rwd_Inds = [];
noRwd_Inds = [];
for i = 1:length(s.behSessionData)
    if s.behSessionData(i).rewardR == 1 || s.behSessionData(i).rewardL == 1
        rwd_Inds = [rwd_Inds i];
    end
    if s.behSessionData(i).rewardR == 0 || s.behSessionData(i).rewardL == 0
        noRwd_Inds = [noRwd_Inds i];
    end
end

%set max trial number for plotting
yLimMax = max([length(blockR) length(blockL) length(lickR_Inds) length(lickL_Inds)]);

rwd_Inds = s.responseInds(ismember(s.responseInds, rwd_Inds));
noRwd_Inds = s.responseInds(~ismember(s.responseInds, rwd_Inds));

%lick side x outcome
[~,lickRrwd_Inds,~] = intersect(rwd_Inds, lickR_Inds);
[~,lickRnoRwd_Inds,~] = intersect(noRwd_Inds, lickR_Inds);
[~,lickLrwd_Inds,~] = intersect(rwd_Inds, lickL_Inds);
[~,lickLnoRwd_Inds,~] = intersect(noRwd_Inds, lickL_Inds);

changeChoice = [false abs(diff(s.allChoices)) > 0];
changeChoice_Inds = find(changeChoice == 1);
stayChoice_Inds = find(changeChoice == 0);


%% smooth rewards over time
sessionTime = [sessionData(1).CSon:sessionData(end).CSon + 3000] - sessionData(1).CSon;       %pad time for reward on last trial

sessionRwds = [sessionData.rewardTime] - sessionData(1).CSon;     %baseline to start time and convert to s from ms
session_rwdsArray = zeros(1,length(sessionTime));
sessionRwds = sessionRwds(rwd_Inds);
session_rwdsArray(sessionRwds) = 1;

boxKern = ones(1,60000);                                       %smooth rewards over time
sessionRwdsSmooth = conv(session_rwdsArray, boxKern);
sessionRwdsSmooth = sessionRwdsSmooth(1:(end-(length(boxKern)-1)));


%% Generate smoothed choice-history values       
[choiceHx_Sorted,choiceHx_Inds] = sort(s.choiceHx);                                          %sort by choice history
choiceHx_Lims = [find(choiceHx_Sorted > -0.9999, 1) find(choiceHx_Sorted > 0.9999, 1)];           %find where choice hist is all L or all R


%% Generate smoothed reward-history values over trials

[rwdHx_Sorted,tmpInds] = sort(s.rwdHx);
rwdHx_Inds = s.responseInds(tmpInds);

%outcome indices for rwd hist
[~,rwdHxRwd_Inds,~] = intersect(rwdHx_Inds, rwd_Inds); 
[~,rwdHxNoRwd_Inds,~] = intersect(rwdHx_Inds, noRwd_Inds); 

%for tercile analysis
tercile = floor(length(rwdHx_Inds)/3);
rwdHxI_Inds = rwdHx_Inds(1:tercile);
rwdHxII_Inds = rwdHx_Inds(tercile+1:tercile*2);
rwdHxIII_Inds = rwdHx_Inds(tercile*2+1:end);


%outcome indices for rwd hist divisions
[~,rwdHxIrwd_Inds,~] = intersect(rwdHxI_Inds, rwd_Inds);          [~,rwdHxInoRwd_Inds,~] = intersect(rwdHxI_Inds, noRwd_Inds);
[~,rwdHxIIrwd_Inds,~] = intersect(rwdHxII_Inds, rwd_Inds);          [~,rwdHxIInoRwd_Inds,~] = intersect(rwdHxII_Inds, noRwd_Inds);
[~,rwdHxIIIrwd_Inds,~] = intersect(rwdHxIII_Inds, rwd_Inds);          [~,rwdHxIIInoRwd_Inds,~] = intersect(rwdHxIII_Inds, noRwd_Inds);

%choice indices for rwd hist divisions
% [~,rwdHxIR_Inds,~] = intersect(rwdHxI_Inds, lickR_Inds);          [~,rwdHxIL_Inds,~] = intersect(rwdHxI_Inds, lickL_Inds);
% [~,rwdHxIIR_Inds,~] = intersect(rwdHxII_Inds, lickR_Inds);          [~,rwdHxIIL_Inds,~] = intersect(rwdHxII_Inds, lickL_Inds);
% [~,rwdHxIIIR_Inds,~] = intersect(rwdHxIII_Inds, lickR_Inds);          [~,rwdHxIIIL_Inds,~] = intersect(rwdHxIII_Inds, lickL_Inds);


%% generate indeces for Q learning model values

if p.Results.modelsFlag
    Qmdls = qLearning_fit2LR([sessionName '.asc']);
    pe = [Qmdls.fourParams_twoLearnRates_tForget.pe];
    peSD = movstd(pe, 10);
    %rBar = [Qmdls.fiveParams_opponency.rBar];

    [pe_Sorted,pe_Inds] = sort(pe);
    [peSD_Sorted,peSD_Inds] = sort(peSD(1:end-1));

    peSD_Inds = peSD_Inds + 1;
    pe_Lim = find(pe_Sorted > 0, 1);
    
    
    [rBar_sorted rBar_Inds] = sort(rBar);
    tercile = floor(length(rBar_Inds)/3);
    rBarI_Inds = rBar_Inds(1:tercile);
    rBarII_Inds = rBar_Inds(tercile+1:tercile*2);
    rBarIII_Inds = rBar_Inds(tercile*2+1:end);
end





%% sort licking for sdfs

for j = 1:length(sessionData)
    trialDurDiff(j) = (sessionData(j).trialEnd - sessionData(j).CSon)- p.Results.tf*1000;
end
trialDurDiff(end) = 0;  %to account for no trialEnd timestamp on last trial

%initialize lick matrices
allTrial_lickMatx = NaN(length(sessionData),length(time)); 
rwd_lickMatx = NaN(length(rwd_Inds), length(time));
noRwd_lickMatx = NaN(length(noRwd_Inds), length(time));
choiceR_lickMatx = NaN(length(lickR_Inds),length(time));
choiceL_lickMatx = NaN(length(lickL_Inds),length(time));
rwdHxI_lickMatx = NaN(length(rwdHxI_Inds), length(time));
rwdHxII_lickMatx = NaN(length(rwdHxII_Inds), length(time));
rwdHxIII_lickMatx = NaN(length(rwdHxIII_Inds), length(time));


for j = 1:length(allTrial_lick)
    tempLick = allTrial_lick{1,j};
    tempLick = tempLick + p.Results.tb*1000; % add this to pad time for SDF
    allTrial_lickMatx(j,tempLick) = 1;
    if trialDurDiff(j) < 0
        allTrial_lickMatx(j, isnan(allTrial_lickMatx(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
    else
        allTrial_lickMatx(j, isnan(allTrial_lickMatx(j,:))) = 0;
    end
    if sum(allTrial_lickMatx(j,:)) == 0     %if there is no spike data for this trial, don't count it
        allTrial_lickMatx(j,:) = NaN;
    end
end

%sort by receipt of reward
for j = 1:length(rwd_Inds)
    rwd_lickMatx(j,:) = allTrial_lickMatx(rwd_Inds(j),:);
end
for j = 1:length(noRwd_Inds)
    noRwd_lickMatx(j,:) = allTrial_lickMatx(noRwd_Inds(j),:);
end

%sort by lick choice
for j = 1:length(lickR_Inds)
    choiceR_lickMatx(j,:) = allTrial_lickMatx(lickR_Inds(j),:);
end
for j = 1:length(lickL_Inds)
    choiceL_lickMatx(j,:) = allTrial_lickMatx(lickL_Inds(j),:);
end

%sort licks by reward history
for j = 1:length(rwdHxI_Inds)
    rwdHxI_lickMatx(j,:) = allTrial_lickMatx(rwdHxI_Inds(j),:);
end
for j = 1:length(rwdHxII_Inds)
    rwdHxII_lickMatx(j,:) = allTrial_lickMatx(rwdHxII_Inds(j),:);
end
for j = 1:length(rwdHxIII_Inds)
    rwdHxIII_lickMatx(j,:) = allTrial_lickMatx(rwdHxIII_Inds(j),:);
end


%% sort spikes for trial averages
for i = 1:length(clust)
    % Initialize matrices for SDF
    
    allTrial_spikeMatx = NaN(length(sessionData),length(time));         
    CSplus_Index = find(CSplus_Inds);                                   
    CSplus_spikeMatx = NaN(length(CSplus_Index),length(time));
    CSminus_Index = find(CSminus_Inds);
    CSminus_spikeMatx = NaN(length(CSminus_Index),length(time)); 
    choiceR_spikeMatx = NaN(length(lickR_Inds),length(time));
    choiceL_spikeMatx = NaN(length(lickL_Inds),length(time)); 
    rwd_spikeMatx = NaN(length(rwd_Inds), length(time));
    noRwd_spikeMatx = NaN(length(noRwd_Inds), length(time));
    rwdHxI_spikeMatx = NaN(length(rwdHxI_Inds), length(time));
    rwdHxII_spikeMatx = NaN(length(rwdHxII_Inds), length(time));
    rwdHxIII_spikeMatx = NaN(length(rwdHxIII_Inds), length(time));
%     rBarI_spikeMatx = NaN(length(rBarI_Inds), length(time));
%     rBarII_spikeMatx = NaN(length(rBarII_Inds), length(time));
%     rBarIII_spikeMatx = NaN(length(rBarIII_Inds), length(time));
    rwdPrev_spikeMatx = NaN(length(rwd_Inds) - 1,length(time));
    noRwdPrev_spikeMatx = NaN(length(noRwd_Inds) - 1,length(time));
    
    for j = 1:length(allTrial_spike)
        tempSpike = allTrial_spike{i,j};
        tempSpike = tempSpike + p.Results.tb*1000; % add this to pad time for SDF
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
    for j = 1:length(CSplus_Index)
        CSplus_spikeMatx(j,:) = allTrial_spikeMatx(CSplus_Index(j),:);
    end
    for j = 1:length(CSminus_Index)
        CSminus_spikeMatx(j,:) = allTrial_spikeMatx(CSminus_Index(j),:);
    end
    for j = 1:length(lickR_Inds)
        choiceR_spikeMatx(j,:) = allTrial_spikeMatx(lickR_Inds(j),:);
    end
    for j = 1:length(lickL_Inds)
        choiceL_spikeMatx(j,:) = allTrial_spikeMatx(lickL_Inds(j),:);
    end  
    for j = 1:length(rwd_Inds)
        rwd_spikeMatx(j,:) = allTrial_spikeMatx(rwd_Inds(j),:);
    end
    for j = 1:length(noRwd_Inds)
        noRwd_spikeMatx(j,:) = allTrial_spikeMatx(noRwd_Inds(j),:);
    end
    for j = 1:length(rwd_Inds)-1
        rwdPrev_spikeMatx(j,:) = allTrial_spikeMatx(rwd_Inds(j)+1,:);
    end
    for j = 1:length(noRwd_Inds)-1
        noRwdPrev_spikeMatx(j,:) = allTrial_spikeMatx(noRwd_Inds(j)+1,:);
    end
    
    choiceTrial_spikeMatx = allTrial_spikeMatx(s.responseInds,:);
    for j = 1:length(changeChoice_Inds)
        changeChoice_spikeMatx(j,:) = choiceTrial_spikeMatx(changeChoice_Inds(j),:);
    end
    for j = 1:length(stayChoice_Inds)
        stayChoice_spikeMatx(j,:) = choiceTrial_spikeMatx(stayChoice_Inds(j),:);
    end   
        
    %initialize matrices for sorting by rwd hist and outcome
    rwdHxIRwd_spikeMatx = [];       rwdHxInoRwd_spikeMatx = [];         rwdHxIR_spikeMatx = [];      rwdHxIL_spikeMatx = [];
    rwdHxIIRwd_spikeMatx = [];      rwdHxIInoRwd_spikeMatx = [];        rwdHxIIR_spikeMatx = [];     rwdHxIIL_spikeMatx = [];
    rwdHxIIIRwd_spikeMatx = [];     rwdHxIIInoRwd_spikeMatx = [];       rwdHxIIIR_spikeMatx = [];    rwdHxIIIL_spikeMatx = [];
    
   
    %sort spikes by reward history
    for j = 1:length(rwdHxI_Inds)
        rwdHxI_spikeMatx(j,:) = allTrial_spikeMatx(rwdHxI_Inds(j),:);
    end
    for j = 1:length(rwdHxII_Inds)
        rwdHxII_spikeMatx(j,:) = allTrial_spikeMatx(rwdHxII_Inds(j),:);
    end
    for j = 1:length(rwdHxIII_Inds)
        rwdHxIII_spikeMatx(j,:) = allTrial_spikeMatx(rwdHxIII_Inds(j),:);
    end

%     for j = 1:length(rBarI_Inds)
%         rBarI_spikeMatx(j,:) = allTrial_spikeMatx(rBarI_Inds(j),:);
%     end
%     for j = 1:length(rBarII_Inds)
%         rBarII_spikeMatx(j,:) = allTrial_spikeMatx(rBarII_Inds(j),:);
%     end
%     for j = 1:length(rBarIII_Inds)
%         rBarIII_spikeMatx(j,:) = allTrial_spikeMatx(rBarIII_Inds(j),:);
%     end
    

    % Downsample raster if necessary
    tempSpike = [sessionData.(spikeFields{clust(i)})];
    avgFiringRate = length(tempSpike)/((tempSpike(end)-tempSpike(1))/1000);

    downsampFlag = false;
    if avgFiringRate >= 12
        downsampFactor = floor(avgFiringRate/6);
        downsampFlag = true;
        for j = 1:size(allTrial_spike,2)
            allTrial_spike{i,j} = downsample(allTrial_spike{i,j}, downsampFactor);
        end   
    end
 
    
    %% lick latency and recent rwd hist comparison
    
    smoothWin = 250;
    trialBeg = p.Results.tb*1000;
    trialEnd = p.Results.tb*1000 + 1500;

    for j = 1:length(s.responseInds)
        spikeTemp = fastsmooth(allTrial_spikeMatx(s.responseInds(j),:)*1000, smoothWin, 3);
        maxFRtimeTemp = find(spikeTemp == max(spikeTemp(trialBeg:trialEnd)), 1);
        preTrialSpikeCount(1,j) = sum(allTrial_spikeMatx(s.responseInds(j), 1:trialBeg));
        if isempty(maxFRtimeTemp)
            maxFRtime(1,j) = NaN;
        else
            maxFRtime(1,j) = maxFRtimeTemp;
        end
        maxFRtrial(1,j) = max(spikeTemp(trialBeg:trialEnd));
    end
    % zscore fr's and remove bad licks
    maxFRtrialRsp = zscore(maxFRtrial); maxFRtrialRsp = maxFRtrialRsp(s.lickLatInds);
    maxFRtimeRsp = maxFRtime(s.lickLatInds);
    preTrialSpikeCountRsp = zscore(preTrialSpikeCount); preTrialSpikeCountRsp = preTrialSpikeCountRsp(s.lickLatInds);
    
    scatters = figure; hold on
    
    mdl = fitlm(maxFRtrialRsp, s.responseLat);
    rSqr = mdl.Rsquared(1).Ordinary(1);
    subplot(2,3,1)
    mdl.plot
    xlabel('Z-Scored Maximum FR')
    ylabel('Z-Scored Lick Latency')
    plotTmp = gca;
    set(plotTmp.Title, 'String', '');
    plotTmp.Children(4).Marker = 'o'; plotTmp.Children(3).LineWidth = 3;
    legend([plotTmp.Children(3)],sprintf(['R^2: ' num2str(rSqr)]))
    
    mdl = fitlm(preTrialSpikeCountRsp, s.responseLat);                             
    rSqr = mdl.Rsquared(1).Ordinary(1);
    subplot(2,3,2)
    mdl.plot
    xlabel('Z-Scored Pre-CS Spike Count')
    ylabel('Z-Scored Lick Latency')
    plotTmp = gca;
    set(plotTmp.Title, 'String', '');
    plotTmp.Children(4).Marker = 'o'; plotTmp.Children(3).LineWidth = 3;
    legend([plotTmp.Children(3)],sprintf(['R^2: ' num2str(rSqr)]))
    
    mdl = fitlm(maxFRtimeRsp, s.responseLat);                             
    rSqr = mdl.Rsquared(1).Ordinary(1);
    subplot(2,3,3)
    mdl.plot
    xlabel('Time of Maximum FR')
    ylabel('Z-Scored Lick Latency')
    plotTmp = gca;
    set(plotTmp.Title, 'String', '');
    plotTmp.Children(4).Marker = 'o'; plotTmp.Children(3).LineWidth = 3;
    legend([plotTmp.Children(3)],sprintf(['R^2: ' num2str(rSqr)]))
    

    subplot(2,3,4); hold on;
    rwdTimeMatx_R = s.rwdTimeMatx;
    rwdTimeMatx_R(find(rwdTimeMatx_R < 0)) = 0;
    mdl = fitlm([rwdTimeMatx_R]', maxFRtrial);
    relevInds = 2:p.Results.timeBins+1;
    coefVals = mdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:p.Results.timeBins)*s.timeBinSize/1000),coefVals,errorL,errorU,'Color', [1 0 0],'linewidth',2)
    
    rwdTimeMatx_L = s.rwdTimeMatx;
    rwdTimeMatx_L(find(rwdTimeMatx_L > 0)) = 0;
    rwdTimeMatx_L(find(rwdTimeMatx_L < 0)) = rwdTimeMatx_L(find(rwdTimeMatx_L < 0)) * -1;
    mdl = fitlm([rwdTimeMatx_L]', maxFRtrial);
    coefVals = mdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:p.Results.timeBins)*s.timeBinSize/1000),coefVals,errorL,errorU,'Color', [0 0 1],'linewidth',2)
    
    mdl = fitlm([s.rwdTimeMatxBin]', maxFRtrial);
    coefVals = mdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:p.Results.timeBins)*s.timeBinSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlabel('reward n seconds back')
    ylabel('\beta Coefficient')
    title('LRM: rewards in time on max fr')
    xlim([0 (p.Results.timeBins*s.timeBinSize/1000 + 5)])
    legend('R only', 'L only', 'all')
    
    subplot(2,3,5); hold on;
    mdl = fitlm([rwdTimeMatx_R]', preTrialSpikeCount);
    coefVals = mdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:p.Results.timeBins)*s.timeBinSize/1000),coefVals,errorL,errorU,'Color', [1 0 0],'linewidth',2)
    
    mdl = fitlm([rwdTimeMatx_L]', preTrialSpikeCount);
    coefVals = mdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:p.Results.timeBins)*s.timeBinSize/1000),coefVals,errorL,errorU,'Color', [0 0 1],'linewidth',2)
    
    mdl = fitlm([s.rwdTimeMatxBin]', preTrialSpikeCount);
    coefVals = mdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:p.Results.timeBins)*s.timeBinSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlabel('reward n seconds back')
    ylabel('\beta Coefficient')
    title('LRM: rewards in time on pre-CS fr')
    xlim([0 (p.Results.timeBins*s.timeBinSize/1000 + 5)])
    legend('R only', 'L only', 'all')
    
   subplot(2,3,6)
   histogram([s.lickLat(s.responseInds) - (maxFRtime-trialBeg)], 'FaceColor', 'c')
   xlabel('lick latency - time of max fr')
   ylabel('counts')
   
   
    set(scatters, 'Position', get(0,'Screensize'))
    cellName = erase(spikeFields{clust(i)}, '_');
    suptitle([sessionName ': ' cellName])

    

    %% session spike averages
    
    %spikes
    sessionSpikes = sessionData(i).allSpikes - sessionData(1).CSon;
    sessionSpikeInd = sessionSpikes(sessionSpikes > 0);
    session_spikeArray = zeros(1,length(sessionTime));      
    session_spikeArray(sessionSpikeInd) = 1000;
    
    boxKern = ones(1,5000) / 5000;                                      smoothKern = ones(1,30000) / 30000;
    sessionSpikeSDF = conv(session_spikeArray, boxKern);                sessionSpikeSDFsmooth = conv(session_spikeArray, smoothKern);
    sessionSpikeSDF = sessionSpikeSDF(1:(end-length(boxKern)-1));       sessionSpikeSDFsmooth = sessionSpikeSDFsmooth(1:(end-length(boxKern)-1));

    trialKern = ones(1,250) / 250;
    trialSpikeSDF = conv(session_spikeArray, trialKern);
    trialSpikeSDF = trialSpikeSDF(1:(end-length(trialKern)-1));

    if length(sessionSpikeSDF) < length(sessionTime)
        sessionSpikeSDF = [sessionSpikeSDF nan(1,length(sessionTime)-length(sessionSpikeSDF))];
        sessionSpikeSDFsmooth = [sessionSpikeSDFsmooth nan(1,length(sessionTime)-length(sessionSpikeSDFsmooth))];
        trialSpikeSDF = [trialSpikeSDF nan(1,length(sessionTime)-length(trialSpikeSDF))];
    end
   

    
    
    %% Plotting
    
    rasters = figure; hold on
    axisColor = [0 0 0];
    set(rasters,'defaultAxesColorOrder',[axisColor; axisColor]);

    % All trials
    r(1) = subplot(8,7,[3 10 17 24]); t(1) = title('All Trials');
    plotSpikeRaster(allTrial_spike(i,:),'PlotType','vertline'); hold on
    plot(repmat([-5000 10000],length(s.blockSwitch),1)', [s.blockSwitch; s.blockSwitch],'r')
    
    % CS+ trials
    r(2) = subplot(8,7,[2 9 16]); t(2) = title('CS Plus');
    plotSpikeRaster(allTrial_spike(i,CSplus_Inds),'PlotType','vertline'); hold on
    plot(repmat([-5000 10000],length(blockSwitch_CSminCorrected),1)', [blockSwitch_CSminCorrected; blockSwitch_CSminCorrected],'r')
    
    
    % CS- trials
    r(3) = subplot(8,7,[23]); t(3) = title('CS Minus');
    if ~all(cellfun(@isempty,allTrial_spike(i,CSminus_Inds)))
        plotSpikeRaster(allTrial_spike(i,CSminus_Inds),'PlotType','vertline'); hold on
        plot([-5000 5000],[sum(CSminus_Inds) sum(CSminus_Inds)],'Color',[192 192 192]/255)
        ylim([0 round(sum(CSplus_Inds)/3)]);
    end
    
    
    % Block_R trials
    r(4) = subplot(8,7,[4 11]); t(4) = title('higher prob on R');
    if ~all(cellfun(@isempty,allTrial_spike(i,blockR)))
        plotSpikeRaster(allTrial_spike(i,blockR),'PlotType','vertline'); hold on
        plot([-5000 5000],[length(blockR) length(blockR)],'Color',[192 192 192]/255)
        bR_switch = find(diff(blockR) > 1);
        plot(repmat([-5000 10000],length(bR_switch),1)', [bR_switch; bR_switch],'r')
        set(gca,'Xticklabel',[]);
    end
    
    
    % Block_L trials
    r(5) = subplot(8,7,[18 25]); t(5) = title('higher prob on L');
    if ~all(cellfun(@isempty,allTrial_spike(i,blockL)))
        plotSpikeRaster(allTrial_spike(i,blockL),'PlotType','vertline'); hold on
        plot([-5000 5000],[length(blockL) length(blockL)],'Color',[192 192 192]/255)
        bL_switch = find(diff(blockL) > 1);
        plot(repmat([-5000 10000],length(bL_switch),1)', [bL_switch; bL_switch],'r')
        set(gca,'Xticklabel',[]);
    end
    
    % R lick broken by rwd vs no rwd
    r(6) = subplot(8,7,[5 12]); t(6) = title('R lick; no rwd v. rwd');
    tmp = find(Rlicks);
    R_rwd = tmp(ismember(tmp,rwd_Inds));
    R_norwd = tmp(ismember(tmp, noRwd_Inds));
    if ~all(cellfun(@isempty,allTrial_spike(i,[R_norwd R_rwd])))
        plotSpikeRaster(allTrial_spike(i,[R_norwd R_rwd]),'PlotType','vertline'); hold on
        plot([-5000 10000],[length(R_norwd) length(R_norwd)],'r')
        plot([-5000 10000],[length([R_norwd R_rwd]) length([R_norwd R_rwd])],'Color',[192 192 192]/255)
        set(gca,'Xticklabel',[]);
    end
    
    % L lick broken by rwd vs no rwd
    r(7) = subplot(8,7,[19 26]); t(7) = title('L lick; no rwd v. rwd');
    tmp = find(Llicks);
    L_rwd = tmp(ismember(tmp,rwd_Inds));
    L_norwd = tmp(ismember(tmp, noRwd_Inds));
    if ~all(cellfun(@isempty,allTrial_spike(i,[L_norwd L_rwd])))
        plotSpikeRaster(allTrial_spike(i,[L_norwd L_rwd]),'PlotType','vertline'); hold on
        plot([-5000 10000],[length(L_norwd) length(L_norwd)],'r')
        plot([-5000 10000],[length([L_norwd L_rwd]) length([L_norwd L_rwd])],'Color',[192 192 192]/255)
        set(gca,'Xticklabel',[]);
    end
    
    % R lick sorted by latency
    r(8) = subplot(8,7,[6 13]); t(8) = title('R licks');
    if ~all(cellfun(@isempty,allTrial_spike(i,lickR_Inds)))
        plotSpikeRaster(allTrial_spike(i,lickR_Inds),'PlotType','vertline'); hold on
        plot([-5000 10000],[length(lickR_Inds) length(lickR_Inds)],'Color',[192 192 192]/255)
        set(gca,'Xticklabel',[]);
    end
        
    % L lick sorted by latency
    r(9) = subplot(8,7,[20 27]); t(9) = title('L licks');
    if ~all(cellfun(@isempty,allTrial_spike(i,lickL_Inds)))
        plotSpikeRaster(allTrial_spike(i,lickL_Inds),'PlotType','vertline'); hold on  
        plot([-5000 10000],[length(lickL_Inds) length(lickL_Inds)],'Color',[192 192 192]/255)
        set(gca,'Xticklabel',[]);
    end
    
    % R lick trials sorted by latency and aligned to first lick
    r(10) = subplot(8,7,[7 14]); t(10) = title('R licks aligned');
    subLat = s.lickLat(lickR_Inds);
    tempSpikes = allTrial_spike(i,lickR_Inds);
    if ~all(cellfun(@isempty,allTrial_spike(i,lickR_Inds)))
        lickR_spike = arrayfun(@(k) tempSpikes{k}-subLat(k),1:length(subLat),'UniformOutput',false);
        plotSpikeRaster(lickR_spike,'PlotType','vertline'); hold on
        plot([-5000 10000],[length(lickR_Inds) length(lickR_Inds)],'Color',[192 192 192]/255)
        set(gca,'Xticklabel',[]);
    end
    
    % L lick trials sorted by latency and aligned to first lick
    r(11) = subplot(8,7,[21 28]); t(11) = title('L licks aligned');
    subLat = s.lickLat(lickL_Inds);
    tempSpikes = allTrial_spike(i,lickL_Inds);
    if ~all(cellfun(@isempty,allTrial_spike(i,lickL_Inds)))
        lickL_spike = arrayfun(@(k) tempSpikes{k}-subLat(k),1:length(subLat),'UniformOutput',false);
        plotSpikeRaster(lickL_spike,'PlotType','vertline'); hold on
        plot([-5000 10000],[length(lickL_Inds) length(lickL_Inds)],'Color',[192 192 192]/255)
        set(gca,'Xticklabel',[]);
    end
    
    % rwd-history
    r(12) = subplot(8,7,[33 40]); t(12) = title('rwd hist L->H');
    plotSpikeRaster(allTrial_spike(i,rwdHx_Inds),'PlotType','vertline'); hold on
    set(gca,'Xticklabel',[]);

    % rwd-history by outcome
    r(13) = subplot(8,7,[34 41]); t(13) = title('rwd hist L->H; no rwd v. rwd');
    rwdHxOutcome = [allTrial_spike(i,rwdHx_Inds(rwdHxNoRwd_Inds)) allTrial_spike(i,rwdHx_Inds(rwdHxRwd_Inds))];
    plotSpikeRaster(rwdHxOutcome,'PlotType','vertline'); hold on
    plot([-5000 10000],[length(rwdHxNoRwd_Inds) length(rwdHxNoRwd_Inds)],'-r')
    set(gca,'Xticklabel',[]);
     
    % sorted by choice history
    r(14) = subplot(8,7,[35 42]); t(14) = title('choice hist L -> R');
    plotSpikeRaster(allTrial_spike(i,choiceHx_Inds),'PlotType','vertline'); hold on
    plot(repmat([-5000 10000],length(choiceHx_Lims),1)', [choiceHx_Lims; choiceHx_Lims],'r')
    set(gca,'Xticklabel',[]);
    
%     % sorted by prediction error
%     r(15) = subplot(8,7,[48 55]); t(15) = title('pe L -> H');
%     plotSpikeRaster(allTrial_spike(i,pe_Inds),'PlotType','vertline'); hold on
%     plot(repmat([-5000 10000],length(pe_Lim),1)', [pe_Lim; pe_Lim],'r')
%     set(gca,'Xticklabel',[]);
%     
%     % sorted by prediction error variance
%     r(16) = subplot(8,7,[49 56]); t(16) = title('pe sd L -> H');
%     plotSpikeRaster(allTrial_spike(i,peSD_Inds),'PlotType','vertline'); hold on
%     set(gca,'Xticklabel',[]);
    
    % CS+/CS- SDFs
    z(1) = subplot(8,7,[30 37]); hold on
    mySDF_CSplus = fastsmooth(nanmean(CSplus_spikeMatx, 1)*1000, smoothWin);
    mySDF_CSminus = fastsmooth(nanmean(CSminus_spikeMatx, 1)*1000, smoothWin);
    plot(time, mySDF_CSminus(1:length(time)),'-','Color', [0.5 0 0.8],'LineWidth',2); hold on
    plot(time, mySDF_CSplus(1:length(time)),'b','LineWidth',2)
    legend({'CS Minus','CS Plus'},'FontSize',6,'Location','northeast')
    
    % rwd/noRwd and lick rate SDFs
    z(2) = subplot(8,7,[31 38]); hold on
    mySDF_rwdR = fastsmooth(nanmean(rwd_spikeMatx(lickRrwd_Inds,:), 1)*1000, smoothWin);
    mySDF_rwdL = fastsmooth(nanmean(rwd_spikeMatx(lickLrwd_Inds,:), 1)*1000, smoothWin);
    mySDF_noRwdR = fastsmooth(nanmean(noRwd_spikeMatx(lickRnoRwd_Inds,:), 1)*1000, smoothWin);
    mySDF_noRwdL = fastsmooth(nanmean(noRwd_spikeMatx(lickLnoRwd_Inds,:), 1)*1000, smoothWin);
    plot(time, mySDF_rwdR(1:length(time)),'-','Color', [0.5 0 0.8],'LineWidth',1)
    plot(time, mySDF_noRwdR(1:length(time)),'-','Color',[0.9 0.6 1],'LineWidth',1)
    plot(time, mySDF_rwdL(1:length(time)),'-b','LineWidth',1)
    plot(time, mySDF_noRwdL(1:length(time)),'-','Color', [0.6 0.6 1],'LineWidth',1)
    legend({'rwd R','no rwd R','rwd L','no rwd L'},'FontSize',6,'Location','northeast')
    
    % rdwHx terciles SDFs
    z(3) = subplot(8,7,[32 39]); t(5) = title('rwd hist'); hold on
    mySDF_rwdHxIlicks = fastsmooth(nanmean(rwdHxI_lickMatx, 1)*1000, smoothWin);
    mySDF_rwdHxIIlicks = fastsmooth(nanmean(rwdHxII_lickMatx, 1)*1000, smoothWin);
    mySDF_rwdHxIIIlicks = fastsmooth(nanmean(rwdHxIII_lickMatx, 1)*1000, smoothWin);
%    mySDF_rwdHxIVlicks = fastsmooth(nanmean(rwdHxIV_lickMatx, 1)*1000, smoothWin);
    mySDF_rwdHxI = fastsmooth(nanmean(rwdHxI_spikeMatx, 1)*1000, 250);
    mySDF_rwdHxII = fastsmooth(nanmean(rwdHxII_spikeMatx, 1)*1000, 250);
    mySDF_rwdHxIII = fastsmooth(nanmean(rwdHxIII_spikeMatx, 1)*1000, 250);
%    mySDF_rwdHxIV = fastsmooth(nanmean(rwdHxIV_spikeMatx, 1)*1000, 250);
    yyaxis right
%    plot(time, mySDF_rwdHxIVlicks(1:length(time)),'--','Color', [0.6 1 0.6], 'LineWidth',1); hold on 
    plot(time, mySDF_rwdHxIIIlicks(1:length(time)),'-','Color',[0 0 0],'LineWidth',2)
    plot(time, mySDF_rwdHxIIlicks(1:length(time)),'-','Color', [0.4 0.4 0.4],'LineWidth',2)
    plot(time, mySDF_rwdHxIlicks(1:length(time)),'-','Color', [0.8 0.8 0.8],'LineWidth',2)
    ylim([0 (max(mySDF_rwdHxIIIlicks)*5)])
    yyaxis left
%    plot(time, mySDF_rwdHxIV(1:length(time)),'b','LineWidth',2)
    plot(time, mySDF_rwdHxIII(1:length(time)),'b', 'LineWidth',2)
    plot(time, mySDF_rwdHxII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
    plot(time, mySDF_rwdHxI(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
    legend({'High','Middle','Low'},'FontSize',5,'Location','northeast')
    set(gca,'Xticklabel',[]);
    
    %sdf by reward hist, rewarded trials
    z(4) = subplot(8,7,[46 53]); t(5) = title('rwd hist - rwd'); hold on
    mySDF_rwdHxIrwd = fastsmooth(nanmean(rwdHxI_spikeMatx(rwdHxIrwd_Inds,:), 1)*1000, smoothWin);
    mySDF_rwdHxIIrwd = fastsmooth(nanmean(rwdHxII_spikeMatx(rwdHxIIrwd_Inds,:), 1)*1000, smoothWin);
    mySDF_rwdHxIIIrwd = fastsmooth(nanmean(rwdHxIII_spikeMatx(rwdHxIIIrwd_Inds,:), 1)*1000, smoothWin);
%    mySDF_rwdHxIVrwd = fastsmooth(nanmean(rwdHxIV_spikeMatx(rwdHxIVrwd_Inds,:), 1)*1000, 250);
%    plot(time, mySDF_rwdHxIVrwd(1:length(time)),'g','LineWidth',2)
    plot(time, mySDF_rwdHxIIIrwd(1:length(time)),'b', 'LineWidth',2)
    plot(time, mySDF_rwdHxIIrwd(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
    plot(time, mySDF_rwdHxIrwd(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
    
    %sdf by reward hist, NON-rewarded trials
    z(5) = subplot(8,7,[47 54]); t(6) = title('rwd hist - no rwd'); hold on
    mySDF_rwdHxInoRwd = fastsmooth(nanmean(rwdHxI_spikeMatx(rwdHxInoRwd_Inds,:), 1)*1000, smoothWin);
    mySDF_rwdHxIInoRwd = fastsmooth(nanmean(rwdHxII_spikeMatx(rwdHxIInoRwd_Inds,:), 1)*1000, smoothWin);
    mySDF_rwdHxIIInoRwd = fastsmooth(nanmean(rwdHxIII_spikeMatx(rwdHxIIInoRwd_Inds,:), 1)*1000, smoothWin);
%    mySDF_rwdHxIVnoRwd = fastsmooth(nanmean(rwdHxIV_spikeMatx(rwdHxIVnoRwd_Inds,:), 1)*1000, 250);
%    plot(time, mySDF_rwdHxIVnoRwd(1:length(time)),'g','LineWidth',2)
    plot(time, mySDF_rwdHxIIInoRwd(1:length(time)),'b', 'LineWidth',2)
    plot(time, mySDF_rwdHxIInoRwd(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
    plot(time, mySDF_rwdHxInoRwd(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
    
    
%   %sdf by rwd hist, R choice trials
%     z(6) = subplot(8,7,[48 55]); t(5) = title('rwd hist - R choice'); hold on
%     mySDF_rwdHxIR = fastsmooth(nanmean(rwdHxI_spikeMatx(rwdHxIR_Inds,:), 1)*1000, smoothWin);
%     mySDF_rwdHxIIR = fastsmooth(nanmean(rwdHxII_spikeMatx(rwdHxIIR_Inds,:), 1)*1000, smoothWin);
%     mySDF_rwdHxIIIR = fastsmooth(nanmean(rwdHxIII_spikeMatx(rwdHxIIIR_Inds,:), 1)*1000, smoothWin);
%     plot(time, mySDF_rwdHxIIIR(1:length(time)),'b', 'LineWidth',2)
%     plot(time, mySDF_rwdHxIIR(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
%     plot(time, mySDF_rwdHxIR(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
%     
%     %sdf by rwd hist, L choice trials
%     z(7) = subplot(8,7,[49 56]); t(5) = title('rwd hist - L choice'); hold on
%     mySDF_rwdHxIL = fastsmooth(nanmean(rwdHxI_spikeMatx(rwdHxIL_Inds,:), 1)*1000, smoothWin);
%     mySDF_rwdHxIIL = fastsmooth(nanmean(rwdHxII_spikeMatx(rwdHxIIL_Inds,:), 1)*1000, smoothWin);
%     mySDF_rwdHxIIIL = fastsmooth(nanmean(rwdHxIII_spikeMatx(rwdHxIIIL_Inds,:), 1)*1000, smoothWin);
%     plot(time, mySDF_rwdHxIIIL(1:length(time)),'b', 'LineWidth',2)
%     plot(time, mySDF_rwdHxIIL(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
%     plot(time, mySDF_rwdHxIL(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)

    z(8) = subplot(8,7,[48 55]); title('prev trial outcome'); hold on;
    mySDF_rwdPrev = fastsmooth(nanmean(rwdPrev_spikeMatx, 1)*1000, smoothWin);
    mySDF_noRwdPrev = fastsmooth(nanmean(noRwdPrev_spikeMatx, 1)*1000, smoothWin);
    plot(time, mySDF_rwdPrev(1:length(time)),'-b','LineWidth',2)
    plot(time, mySDF_noRwdPrev(1:length(time)),'-','Color',[0.7 0 1],'LineWidth',2)
    legend('rwd', 'no rwd')
    
    z(9) = subplot(8,7,[49 56]); hold on;
    mySDF_stay = fastsmooth(nanmean(stayChoice_spikeMatx)*1000, smoothWin);
    mySDF_change = fastsmooth(nanmean(changeChoice_spikeMatx, 1)*1000, smoothWin);
    plot(time, mySDF_stay(1:length(time)),'-b','LineWidth',2)
    plot(time, mySDF_change(1:length(time)),'-','Color',[0.7 0 1],'LineWidth',2)
    legend('stay', 'switch')
    
    %ISI histogram
    isi = subplot(8,7,[29 36]); t(15) = title('ISI Histogram'); hold on
    histogram(diff([sessionData.(spikeFields{clust(i)})]),0:100);
    xlim([0 20])
    
    %reward and choice behavior
    slopeSpace = subplot(8,7,[1 8 15]); t(16) = title('Behavior'); hold on
    set(slopeSpace,'YDir','reverse')
    plot(s.choiceSlope,1:length(s.choiceSlope),'b','linewidth',1.5)
    plot(s.rwdSlope,1:length(s.rwdSlope),'k','linewidth',1.5)
    xlim([0 90])
    ylim([0 sum(CSplus_Inds)])
    xlabel('<-- Left Choice (Slope) Right Choice -->')
    ylabel('Trials')
    legend({'C','R'},'FontSize',5,'location','best')
    plot(repmat([-5000 10000],length(blockSwitch_CSminCorrected),1)', [blockSwitch_CSminCorrected; blockSwitch_CSminCorrected],'r')
    
    %rwd hist and firing rate across session
    subplot(8,7,[43 44 45; 50 51 52])
    yyaxis left
    plot(sessionSpikeSDF,'Color', [0.5 0.5 0.5], 'LineWidth',2); hold on;
    plot(sessionSpikeSDFsmooth,'-','Color', [0.5 0 0.8], 'LineWidth',2); hold on;
    ylabel('Spikes/s')
    xlim([60000 length(sessionTime)])
    yyaxis right
    plot(sessionRwdsSmooth, 'b', 'LineWidth', 2);
    legend('fr', 'fr smooth', 'rwd hist', 'location','best')
    xlabel('Time')
    set(gca,'Xtick',[]);
    xlim([60000 length(sessionTime)])
    rMin = ylim;
    rMag = rMin(2)*0.9;
    for j = 1:length(s.behSessionData)
        if s.behSessionData(j).rewardL == 1 | s.behSessionData(j).rewardR == 1
            xTemp = s.behSessionData(j).rewardTime - s.behSessionData(1).CSon;
            plot([xTemp xTemp],[rMag rMin(2)],'b')
        end
    end
    
    %determine minimum and maximum y values for SDFs
    minY = min([mySDF_rwdHxI(trialBeg:trialEnd) mySDF_rwdHxII(trialBeg:trialEnd) mySDF_rwdHxIII(trialBeg:trialEnd)]);
    maxY = max([mySDF_rwdHxI(trialBeg:trialEnd) mySDF_rwdHxII(trialBeg:trialEnd) mySDF_rwdHxIII(trialBeg:trialEnd)]);
    if minY - 0.07*minY > 0
        SDFyLimMin = minY - 0.07*minY;
    else
        SDFyLimMin = 0;
    end
    SDFyLimMax = maxY + 0.07*maxY;
    
    linkaxes([r z],'x')
    xlim(r(1),[min(time)+500 max(time)-500])
    ylim(r(5),[0 yLimMax])
    linkaxes(r(5:14),'y')
    ylim(z(5),[SDFyLimMin SDFyLimMax])
    linkaxes(z(1:9),'y')

    
    % max size
    set(rasters, 'Position', get(0,'Screensize'))
    
    
    axes( 'Position', [0, 0.95, 1, 0.05] ) ; % set axes for the 'text' call below
    set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
    if downsampFlag == true
        text( 0.5, 0, [sessionName ': ' spikeFields{clust(i)} '. Downsampled By ' num2str(downsampFactor)], 'FontSize', 14', 'FontWeight', 'Bold', ...
          'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom','Interpreter','none') ;
    else
        text( 0.5, 0, [sessionName ': ' spikeFields{clust(i)}], 'FontSize', 14', 'FontWeight', 'Bold', ...
          'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom','Interpreter','none') ;
    end
    
    
    if p.Results.saveFigFlag == 1 
        saveFigurePDF(rasters,[savepath sep sessionName '_' spikeFields{clust(i)}])
        saveFigurePDF(scatters,[savepath sep sessionName '_' spikeFields{clust(i)} '_scatter'])
    end
end