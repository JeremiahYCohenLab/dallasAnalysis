function [mdlStruct] = spikeAnalysis_opMD(sessionName, varargin)


p = inputParser;
% default parameters if none given
p.addParameter('intanFlag', 0);
p.addParameter('revForFlag', 0);
p.addParameter('modelsFlag', 1);
p.addParameter('biasFlag',0);
p.addParameter('cellName', ['all']);
p.addParameter('trialList', []);
p.parse(varargin{:});

cellName = p.Results.cellName;

% Path
[root,sep] = currComputer();

[animalName] = strtok(sessionName, 'd');
animalName = animalName(2:end);

if isstrprop(sessionName(end), 'alpha')
    sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
else
    sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
end
sortedFolder = dir(sortedFolderLocation);


if p.Results.intanFlag
    if any(~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'))) == 1
        sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'));
        load([sortedFolderLocation sortedFolder(sessionDataInd).name])
    else
        [sessionData] = generateSessionData_intan_operantMatching(sessionName);
    end
else
    if any(~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'))) == 1
        sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'));
        load([sortedFolderLocation sortedFolder(sessionDataInd).name])
    else
        [sessionData] = generateSessionData_nL_operantMatching(sessionName);
    end
end
[s] = behAnalysisNoPlot_opMD(sessionName, 'revForFlag', p.Results.revForFlag);


if p.Results.modelsFlag
    Qmdls = qLearning_fitOpponency([sessionName '.asc'], p.Results.revForFlag);
    rBar = Qmdls.fiveParams_opponency_zero.rBar;
%    R = Qmdls.fourParams_gradient.R;
    if p.Results.biasFlag
        sumQ = Qmdls.fiveParams_twoLearnRates_alphaForget_bias.Q(:,1) +  Qmdls.fiveParams_twoLearnRates_alphaForget_bias.Q(:,2);
        diffQ = Qmdls.fiveParams_twoLearnRates_alphaForget_bias.Q(:,1) -  Qmdls.fiveParams_twoLearnRates_alphaForget_bias.Q(:,2);
        pe = Qmdls.fiveParams_twoLearnRates_alphaForget_bias.pe;
    else
        sumQ = Qmdls.fourParams_twoLearnRates_alphaForget.Q(:,1) +  Qmdls.fourParams_twoLearnRates_alphaForget.Q(:,2);
        diffQ = Qmdls.fourParams_twoLearnRates_alphaForget.Q(:,1) -  Qmdls.fourParams_twoLearnRates_alphaForget.Q(:,2);
        pe = Qmdls.fourParams_twoLearnRates_alphaForget.pe;
    end
    confQ = abs(diffQ);
    choiceQ = NaN(length(s.allChoices), 1);
    choiceQ(s.allChoices == 1) = Qmdls.fourParams_twoLearnRates_alphaForget.Q((s.allChoices == 1),1);
    choiceQ(s.allChoices == -1) = Qmdls.fourParams_twoLearnRates_alphaForget.Q((s.allChoices == -1),2);
    pe = [0 pe(1:end-1)];
end

if ~isempty(p.Results.trialList)
    s.responseInds = s.responseInds(p.Results.trialList);
    s.sessionRwds = s.sessionRwds(p.Results.trialList);
    rwdHx = rwdHx(p.Results.trialList);
    if p.Results.modelsFlag
        rBar = rBar(p.Results.trialList);
%        R = R(p.Results.trialList);
        sumQ = sumQ(p.Results.trialList);
        diffQ = diffQ(p.Results.trialList);
        confQ = confQ(p.Results.trialList);
        choiceQ = choiceQ(p.Results.trialList);
        pe = pe(p.Results.trialList);
    end
end


%% Sort all spikes into a raster-able matrix

%set time window for spike analyses
tb = 1.5;
tf = 5;
time = -1000*tb:1000*tf;

spikeFields = fields(sessionData);
if iscell(cellName)
    for i = 1:length(cellName)
        cellInd(i) = find(~cellfun(@isempty,strfind(spikeFields,cellName{i})));
    end
elseif regexp(cellName, 'all')
    if p.Results.intanFlag
        cellInd = find(~cellfun(@isempty,strfind(spikeFields,'C_')));
    else
        cellInd = find(~cellfun(@isempty,strfind(spikeFields,'SS')) & ~cellfun(@isempty,strfind(spikeFields,'TT')));
    end
else
    cellInd = find(~cellfun(@isempty,strfind(spikeFields,cellName)));
end
   
allTrial_spike = {};
for k = 1:length(sessionData)
    for i = 1:length(cellInd)
        if k == 1
            prevTrial_spike = [];
            currTrial_lick = [];
        else
            prevTrial_spikeInd = sessionData(k-1).(spikeFields{cellInd(i)}) > sessionData(k-1).trialEnd-tb*1000;
            prevTrial_spike = sessionData(k-1).(spikeFields{cellInd(i)})(prevTrial_spikeInd) - sessionData(k).CSon;
        end
        
        currTrial_spikeInd = sessionData(k).(spikeFields{cellInd(i)}) < sessionData(k).CSon+tf*1000;
        currTrial_spike = sessionData(k).(spikeFields{cellInd(i)})(currTrial_spikeInd) - sessionData(k).CSon;
        
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


%% set time window and smoothing parameters, run analysis for all cells

smoothWin = 250;
trialBeg = tb*1000;
CSoff = tb*1000 + 500;
for i = 1:length(s.behSessionData)
    trialDurDiff(i) = (s.behSessionData(i).trialEnd - s.behSessionData(i).CSon)- tf*1000;
end

mdlStruct = struct;

for i = 1:length(cellInd)
    allTrial_spikeMatx = NaN(length(sessionData),length(time));
    for j = 1:length(allTrial_spike)
        tempSpike = allTrial_spike{i,j};
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
    
%% plot effect of previous reward
    prevRwd_spikeMatx = [];
    prevNoRwd_spikeMatx = [];

    for j = 2:length(s.responseInds)
        if s.allRewardsBinary(j-1) == 1
            prevRwd_spikeMatx = [prevRwd_spikeMatx; allTrial_spikeMatx(s.responseInds(j),:)];
        else
            prevNoRwd_spikeMatx = [prevNoRwd_spikeMatx; allTrial_spikeMatx(s.responseInds(j),:)];
        end
    end

    mySDF_prevRwd = fastsmooth(nanmean(prevRwd_spikeMatx, 1)*1000, smoothWin)';
    mySDF_prevNoRwd = fastsmooth(nanmean(prevNoRwd_spikeMatx, 1)*1000, smoothWin)';
    mySD_prevRwd = fastsmooth(nanstd(prevRwd_spikeMatx, 1)*1000 / sqrt(size(prevRwd_spikeMatx,2)), smoothWin)';
    mySD_prevNoRwd = fastsmooth(nanstd(prevNoRwd_spikeMatx, 1)*1000 / sqrt(size(prevRwd_spikeMatx,2)), smoothWin)';
   
    a = [[-500:2000]' mySDF_prevRwd(1000:3500) mySD_prevRwd(1000:3500) mySDF_prevNoRwd(1000:3500) mySD_prevNoRwd(1000:3500)];
    for j = 1:5
        b(:,j) = downsample(a(:,j),5);
    end
    figure; hold on;
    plot(mySDF_prevNoRwd,'-','Color', [0.5 0 0.8],'LineWidth',2); hold on
    plot(mySDF_prevRwd,'b','LineWidth',2)
    legend('prev no rwd', 'prev rwd');
    
    noResp_spikeMatx = [];
    noRespInds = find(isnan([s.behSessionData.rewardTime]));
    
    if ~isempty(noRespInds)
        for j = 1:length(noRespInds)
            noResp_spikeMatx = [noResp_spikeMatx allTrial_spikeMatx(noRespInds(j),:)];
        end
        mySDF_noResp = fastsmooth(nanmean(noResp_spikeMatx, 1)*1000, smoothWin)';
        figure; plot(mySDF_noResp, '-b', 'linewidth', 2)
        xlim([0 (tb*1000+2000)]);
    end
    
   
    
%% find features of spike rate on each trial    
    for j = 1:length(allTrial_spike)
        if ~isempty(allTrial_spikeMatx(i,j))
            preCSspikeCount(i,j) = sum(allTrial_spikeMatx(j, 1:trialBeg));              %find total spikes before CS on
            postCSspikeCount(i,j) = sum(allTrial_spikeMatx(j, trialBeg:CSoff));

            spikeTemp = fastsmooth(allTrial_spikeMatx(j,:)*1000, smoothWin, 3);         %smooth raw spikes to find features of spike rate
            maxFRcs(i,j) = max(spikeTemp(trialBeg:CSoff));
            minFRcs(i,j) = min(spikeTemp(trialBeg:CSoff));
            if ~isnan(maxFRcs(i,j))
                maxFRtime(i,j) = find(spikeTemp == max(spikeTemp(trialBeg:CSoff)), 1);
            else
                maxFRtime(i,j) = NaN;
            end

        else
            preCSspikeCount(i,j) = NaN;
            postCSspikeCount(i,j) = NaN;
            maxFRcs(i,j) = NaN;
            minFRcs(i,j) = NaN;
            maxFRtime(i,j) = NaN;
        end
    end
    
    
    %organize model parameters for comparison to spike data
%     peMdl = fitlm(peSDm', preCSspikeCount(s.responseInds));
    
%     rBar = [Qmdls.fiveParams_opponency.rBar];
    

%% smooth rewards with different kernels and run linear regression models
    %smooth rewards to compare to spike data
    decayConstList = [1 2 4 8 16 32 64 128 256];
    winEnd = decayConstList*10;
    smoothRwds = [];
    tmpFilt = [];
    for j = 1:length(decayConstList)
        tmpFilt(j,:) = exp(-1/decayConstList(j)*(0:winEnd(end)));
        tmpFilt(j,:) = tmpFilt(j,:)/max(tmpFilt(j,:));
        tempSmoothRwds = conv(s.sessionRwds, tmpFilt(j,:));
        smoothRwds(:,j) = tempSmoothRwds(s.choiceTimes);
    end
    if ~isempty(p.Results.trialList)
        smoothRwds = smoothRwds(p.Results.trialList,:);
    end
    
    if p.Results.modelsFlag
        decayConstList = [0.5 1 2 4 8 16];
        winEnd = decayConstList*10;
        smoothPEs = [];
        tmpFilt = [];
        for j = 1:length(decayConstList)
            tmpFilt(j,:) = exp(-1/decayConstList(j)*(0:winEnd(end)));
            tmpFilt(j,:) = tmpFilt(j,:)/max(tmpFilt(j,:));
            tempSmoothPEs = conv(pe, tmpFilt(j,:));
            smoothPEs(:,j) = tempSmoothPEs(1:end-(length(tmpFilt(j,:))-1)); 
        end
        if ~isempty(p.Results.trialList)
            smoothPEs = smoothPEs(p.Results.trialList,:);
        end
    end
    
    
    
    %create models and put them into the structure
    
    if length(cellInd) > 1
        cellNameTemp = cellName{i};
    else
        cellNameTemp = cellName;
    end
    mdlStruct.(cellNameTemp).maxFRtrialRwdHist = fitlm(smoothRwds, maxFRcs(i,s.responseInds));
    mdlStruct.(cellNameTemp).preCSspikeCountRwdHist = fitlm(smoothRwds, preCSspikeCount(i,s.responseInds));
    mdlStruct.(cellNameTemp).postCSspikeCountRwdHist = fitlm(smoothRwds, postCSspikeCount(i,s.responseInds));
    if p.Results.modelsFlag
        mdlStruct.(cellNameTemp).maxFRtrialQ = fitlm([sumQ diffQ confQ choiceQ s.rwdHx' rBar], maxFRcs(i,s.responseInds));
        mdlStruct.(cellNameTemp).postCSspikeCountQ = fitlm([sumQ diffQ confQ choiceQ s.rwdHx' rBar], postCSspikeCount(i,s.responseInds));
        mdlStruct.(cellNameTemp).preCSspikeCountQ = fitlm([sumQ diffQ confQ choiceQ s.rwdHx' rBar], preCSspikeCount(i,s.responseInds));
        [phasicMinP, phasicMin] = min(mdlStruct.(cellNameTemp).maxFRtrialRwdHist.Coefficients.pValue(2:end, 1));
        mdlStruct.(cellNameTemp).rBar = fitlm([smoothRwds(:,phasicMin) rBar], maxFRcs(i,s.responseInds));
        mdlStruct.(cellNameTemp).PEmaxFR = fitlm(smoothPEs, maxFRcs(i,s.responseInds));
        mdlStruct.(cellNameTemp).PEpreCS = fitlm(smoothPEs, preCSspikeCount(i,s.responseInds));
    end
    

    
end

save([sortedFolderLocation '_spikeMdls.mat'], 'mdlStruct');

% if p.Results.modelsFlag
%     [phasicMinP, phasicMin] = min(mdlStruct.(cellNameTemp).postCSspikeCountRwdHist.Coefficients.pValue(2:end, 1));
%     [tonicMinP, tonicMin] = min(mdlStruct.(cellNameTemp).preCSspikeCountRwdHist.Coefficients.pValue(2:end, 1)); 
% 
%     figure; hold on;
%     title([sessionName ' ' cellName ' best fit reward history'])
%     plot(tmpFilt{phasicMin}, '-', 'Color', [0.7 0 1], 'linewidth', 2)
%     plot(tmpFilt{tonicMin}, '-b',  'linewidth', 2)
%     legend(['phasic pVal = ' num2str(phasicMinP)], ['tonic pVal = ' num2str(tonicMinP)])
%     xlabel('time (s)')
% end


    