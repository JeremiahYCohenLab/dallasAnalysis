function [mdlStruct] = spikeQanalysis_opMD(sessionName, varargin)


p = inputParser;
% default parameters if none given
p.addParameter('saveFigFlag', 1);
p.addParameter('intanFlag', 0);
p.addParameter('revForFlag', 0);
p.addParameter('biasFlag',0);
p.addParameter('cellName', ['all']);
p.addParameter('trialList', []);
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
    savePath = [root animalName sep sessionName(1:end-1) sep 'figures' sep 'session ' sessionName(end) sep];
else
    sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
    savePath = [root animalName sep sessionName sep 'figures' sep];
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



Qmdls = qLearning_fitOpponency([sessionName '.asc'], p.Results.revForFlag);
rBar = Qmdls.fiveParams_opponency.rBar;
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


if ~isempty(p.Results.trialList)
    s.responseInds = s.responseInds(p.Results.trialList);
    s.sessionRwds = s.sessionRwds(p.Results.trialList);
    s.rwdHx = s.rwdHx(p.Results.trialList);
    rBar = rBar(p.Results.trialList);
    sumQ = sumQ(p.Results.trialList);
    diffQ = diffQ(p.Results.trialList);
    confQ = confQ(p.Results.trialList);
    choiceQ = choiceQ(p.Results.trialList);
    pe = pe(p.Results.trialList);
    tercile = floor(length(p.Results.trialList)/3);
else
    tercile = floor(length(s.responseInds)/3);
end


[rwdHx_Sorted, rwdHx_Inds] = sort(s.rwdHx);
rwdHxI_Inds = rwdHx_Inds(1:tercile);
rwdHxII_Inds = rwdHx_Inds(tercile+1:tercile*2);
rwdHxIII_Inds = rwdHx_Inds(tercile*2+1:end);

[rBar_Sorted, rBar_Inds] = sort(rBar);
rBarI_Inds = rBar_Inds(1:tercile);
rBarII_Inds = rBar_Inds(tercile+1:tercile*2);
rBarIII_Inds = rBar_Inds(tercile*2+1:end);

[pe_Sorted, pe_Inds] = sort(pe);
peI_Inds = pe_Inds(1:tercile);
peII_Inds = pe_Inds(tercile+1:tercile*2);
peIII_Inds = pe_Inds(tercile*2+1:end);

[sumQ_Sorted, sumQ_Inds] = sort(sumQ);
sumQI_Inds = sumQ_Inds(1:tercile);
sumQII_Inds = sumQ_Inds(tercile+1:tercile*2);
sumQIII_Inds = sumQ_Inds(tercile*2+1:end);

[diffQ_Sorted, diffQ_Inds] = sort(diffQ);
diffQI_Inds = diffQ_Inds(1:tercile);
diffQII_Inds = diffQ_Inds(tercile+1:tercile*2);
diffQIII_Inds = diffQ_Inds(tercile*2+1:end);

[confQ_Sorted, confQ_Inds] = sort(confQ);
confQI_Inds = confQ_Inds(1:tercile);
confQII_Inds = confQ_Inds(tercile+1:tercile*2);
confQIII_Inds = confQ_Inds(tercile*2+1:end);

[choiceQ_Sorted, choiceQ_Inds] = sort(choiceQ);
choiceQI_Inds = choiceQ_Inds(1:tercile);
choiceQII_Inds = choiceQ_Inds(tercile+1:tercile*2);
choiceQIII_Inds = choiceQ_Inds(tercile*2+1:end);



%% Sort all spikes into a raster-able matrix

%set time window for spike analyses
time = -1000*p.Results.tb:1000*p.Results.tf;

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
        else
            prevTrial_spikeInd = sessionData(k-1).(spikeFields{cellInd(i)}) > sessionData(k-1).trialEnd - p.Results.tb*1000;
            prevTrial_spike = sessionData(k-1).(spikeFields{cellInd(i)})(prevTrial_spikeInd) - sessionData(k).CSon;
        end
        
        currTrial_spikeInd = sessionData(k).(spikeFields{cellInd(i)}) < sessionData(k).CSon + p.Results.tf*1000;
        currTrial_spike = sessionData(k).(spikeFields{cellInd(i)})(currTrial_spikeInd) - sessionData(k).CSon;
        
        allTrial_spike{i,k} = [prevTrial_spike currTrial_spike];
    end
end

% sometimes no spikes are considered 1x0 and sometimes they are []
% plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
allTrial_spike(cellfun(@isempty,allTrial_spike)) = {zeros(1,0)}; 


%% set time window and smoothing parameters, run analysis for all cells

smoothWin = 250;
trialBeg = p.Results.tb*1000;
trialEnd = p.Results.tb*1000 + p.Results.tf*1000;
CSoff = p.Results.tb*1000 + 500;
for i = 1:length(sessionData)
    trialDurDiff(i) = (sessionData(i).trialEnd - sessionData(i).CSon)- p.Results.tf*1000;
end
trialDurDiff(end) = 0;  %to account for no trialEnd timestamp on last trial

mdlStruct = struct;

for i = 1:length(cellInd)
    allTrial_spikeMatx = NaN(length(sessionData),length(time));
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

    rwdHxI_spikeMatx = NaN(length(rwdHxI_Inds), length(time));
    rwdHxII_spikeMatx = NaN(length(rwdHxII_Inds), length(time));
    rwdHxIII_spikeMatx = NaN(length(rwdHxIII_Inds), length(time));
    rBarI_spikeMatx = NaN(length(rBarI_Inds), length(time));
    rBarII_spikeMatx = NaN(length(rBarII_Inds), length(time));
    rBarIII_spikeMatx = NaN(length(rBarIII_Inds), length(time));
    sumQI_spikeMatx = NaN(length(sumQI_Inds), length(time));
    sumQII_spikeMatx = NaN(length(sumQII_Inds), length(time));
    sumQIII_spikeMatx = NaN(length(sumQIII_Inds), length(time));
    diffQI_spikeMatx = NaN(length(diffQI_Inds), length(time));
    diffQII_spikeMatx = NaN(length(diffQII_Inds), length(time));
    diffQIII_spikeMatx = NaN(length(diffQIII_Inds), length(time));
    confQI_spikeMatx = NaN(length(confQI_Inds), length(time));
    confQII_spikeMatx = NaN(length(confQII_Inds), length(time));
    confQIII_spikeMatx = NaN(length(confQIII_Inds), length(time));
    choiceQI_spikeMatx = NaN(length(choiceQI_Inds), length(time));
    choiceQII_spikeMatx = NaN(length(choiceQII_Inds), length(time));
    choiceQIII_spikeMatx = NaN(length(choiceQIII_Inds), length(time));
    peI_spikeMatx = NaN(length(peI_Inds), length(time));
    peII_spikeMatx = NaN(length(peII_Inds), length(time));
    peIII_spikeMatx = NaN(length(peIII_Inds), length(time));
    
    for j = 1:length(rwdHxI_Inds)
        rwdHxI_spikeMatx(j,:) = allTrial_spikeMatx(rwdHxI_Inds(j),:);
    end
    for j = 1:length(rwdHxII_Inds)
        rwdHxII_spikeMatx(j,:) = allTrial_spikeMatx(rwdHxII_Inds(j),:);
    end
    for j = 1:length(rwdHxIII_Inds)
        rwdHxIII_spikeMatx(j,:) = allTrial_spikeMatx(rwdHxIII_Inds(j),:);
    end
    for j = 1:length(rBarI_Inds)
        rBarI_spikeMatx(j,:) = allTrial_spikeMatx(rBarI_Inds(j),:);
    end
    for j = 1:length(rBarII_Inds)
        rBarII_spikeMatx(j,:) = allTrial_spikeMatx(rBarII_Inds(j),:);
    end
    for j = 1:length(rBarIII_Inds)
        rBarIII_spikeMatx(j,:) = allTrial_spikeMatx(rBarIII_Inds(j),:);
    end
    for j = 1:length(sumQI_Inds)
        sumQI_spikeMatx(j,:) = allTrial_spikeMatx(sumQI_Inds(j),:);
    end
    for j = 1:length(sumQII_Inds)
        sumQII_spikeMatx(j,:) = allTrial_spikeMatx(sumQII_Inds(j),:);
    end
    for j = 1:length(sumQIII_Inds)
        sumQIII_spikeMatx(j,:) = allTrial_spikeMatx(sumQIII_Inds(j),:);
    end
    for j = 1:length(diffQI_Inds)
        diffQI_spikeMatx(j,:) = allTrial_spikeMatx(diffQI_Inds(j),:);
    end
    for j = 1:length(diffQII_Inds)
        diffQII_spikeMatx(j,:) = allTrial_spikeMatx(diffQII_Inds(j),:);
    end
    for j = 1:length(diffQIII_Inds)
        diffQIII_spikeMatx(j,:) = allTrial_spikeMatx(diffQIII_Inds(j),:);
    end
    for j = 1:length(confQI_Inds)
        confQI_spikeMatx(j,:) = allTrial_spikeMatx(confQI_Inds(j),:);
    end
    for j = 1:length(confQII_Inds)
        confQII_spikeMatx(j,:) = allTrial_spikeMatx(confQII_Inds(j),:);
    end
    for j = 1:length(confQIII_Inds)
        confQIII_spikeMatx(j,:) = allTrial_spikeMatx(confQIII_Inds(j),:);
    end
    for j = 1:length(choiceQI_Inds)
        choiceQI_spikeMatx(j,:) = allTrial_spikeMatx(choiceQI_Inds(j),:);
    end
    for j = 1:length(choiceQII_Inds)
        choiceQII_spikeMatx(j,:) = allTrial_spikeMatx(choiceQII_Inds(j),:);
    end
    for j = 1:length(choiceQIII_Inds)
        choiceQIII_spikeMatx(j,:) = allTrial_spikeMatx(choiceQIII_Inds(j),:);
    end
    for j = 1:length(peI_Inds)
        peI_spikeMatx(j,:) = allTrial_spikeMatx(peI_Inds(j),:);
    end
    for j = 1:length(peII_Inds)
        peII_spikeMatx(j,:) = allTrial_spikeMatx(peII_Inds(j),:);
    end
    for j = 1:length(peIII_Inds)
        peIII_spikeMatx(j,:) = allTrial_spikeMatx(peIII_Inds(j),:);
    end
    
    %%
    sdfFig = figure;
    
    z(1) = subplot(2,4,1); t(1) = title('rwd hist'); hold on
    mySDF_rwdHxI = fastsmooth(nanmean(rwdHxI_spikeMatx, 1)*1000, 250);
    mySDF_rwdHxII = fastsmooth(nanmean(rwdHxII_spikeMatx, 1)*1000, 250);
    mySDF_rwdHxIII = fastsmooth(nanmean(rwdHxIII_spikeMatx, 1)*1000, 250);
    plot(time, mySDF_rwdHxIII(1:length(time)),'b', 'LineWidth',2)
    plot(time, mySDF_rwdHxII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
    plot(time, mySDF_rwdHxI(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
    legend({'High','Middle','Low'},'FontSize',5,'Location','northeast')
    
    z(2) = subplot(2,4,2); t(2) = title('rBar'); hold on
    mySDF_rBarI = fastsmooth(nanmean(rBarI_spikeMatx, 1)*1000, 250);
    mySDF_rBarII = fastsmooth(nanmean(rBarII_spikeMatx, 1)*1000, 250);
    mySDF_rBarIII = fastsmooth(nanmean(rBarIII_spikeMatx, 1)*1000, 250);
    plot(time, mySDF_rBarIII(1:length(time)),'b', 'LineWidth',2)
    plot(time, mySDF_rBarII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
    plot(time, mySDF_rBarI(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
    
    z(3) = subplot(2,4,3); t(3) = title('sumQ'); hold on
    mySDF_sumQI = fastsmooth(nanmean(sumQI_spikeMatx, 1)*1000, 250);
    mySDF_sumQII = fastsmooth(nanmean(sumQII_spikeMatx, 1)*1000, 250);
    mySDF_sumQIII = fastsmooth(nanmean(sumQIII_spikeMatx, 1)*1000, 250);
    plot(time, mySDF_sumQIII(1:length(time)),'b', 'LineWidth',2)
    plot(time, mySDF_sumQII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
    plot(time, mySDF_sumQI(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
    
    z(4) = subplot(2,4,4); t(4) = title('diffQ'); hold on
    mySDF_diffQI = fastsmooth(nanmean(diffQI_spikeMatx, 1)*1000, 250);
    mySDF_diffQII = fastsmooth(nanmean(diffQII_spikeMatx, 1)*1000, 250);
    mySDF_diffQIII = fastsmooth(nanmean(diffQIII_spikeMatx, 1)*1000, 250);
    plot(time, mySDF_diffQIII(1:length(time)),'b', 'LineWidth',2)
    plot(time, mySDF_diffQII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
    plot(time, mySDF_diffQI(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
    
    z(5) = subplot(2,4,5); t(5) = title('confQ'); hold on
    mySDF_confQI = fastsmooth(nanmean(confQI_spikeMatx, 1)*1000, 250);
    mySDF_confQII = fastsmooth(nanmean(confQII_spikeMatx, 1)*1000, 250);
    mySDF_confQIII = fastsmooth(nanmean(confQIII_spikeMatx, 1)*1000, 250);
    plot(time, mySDF_confQIII(1:length(time)),'b', 'LineWidth',2)
    plot(time, mySDF_confQII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
    plot(time, mySDF_confQI(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
    
    z(6) = subplot(2,4,6); t(6) = title('choiceQ'); hold on
    mySDF_choiceQI = fastsmooth(nanmean(choiceQI_spikeMatx, 1)*1000, 250);
    mySDF_choiceQII = fastsmooth(nanmean(choiceQII_spikeMatx, 1)*1000, 250);
    mySDF_choiceQIII = fastsmooth(nanmean(choiceQIII_spikeMatx, 1)*1000, 250);
    plot(time, mySDF_choiceQIII(1:length(time)),'b', 'LineWidth',2)
    plot(time, mySDF_choiceQII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
    plot(time, mySDF_choiceQI(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
    
    z(7) = subplot(2,4,7); t(7) = title('pe'); hold on
    mySDF_peI = fastsmooth(nanmean(peI_spikeMatx, 1)*1000, 250);
    mySDF_peII = fastsmooth(nanmean(peII_spikeMatx, 1)*1000, 250);
    mySDF_peIII = fastsmooth(nanmean(peIII_spikeMatx, 1)*1000, 250);
    plot(time, mySDF_peIII(1:length(time)),'b', 'LineWidth',2)
    plot(time, mySDF_peII(1:length(time)),'-','Color', [0.4 0.4 1],'LineWidth',2)
    plot(time, mySDF_peI(1:length(time)),'-','Color', [0.8 0.8 1],'LineWidth',2)
    
    
    %determine minimum and maximum y values for SDFs
    minmaxY = minmax([mySDF_rwdHxI(trialBeg:trialEnd-1000) mySDF_rwdHxII(trialBeg:trialEnd-1000) mySDF_rwdHxIII(trialBeg:trialEnd-1000)]);
    if minmaxY(1) - 0.07*minmaxY(1) > 0
        SDFyLimMin = minmaxY(1) - 0.07*minmaxY(1);
    else
        SDFyLimMin = 0;
    end
    SDFyLimMax = minmaxY(2) + 0.07*minmaxY(2);
    
    linkaxes(z,'xy')
    xlim(z(1),[min(time)+500 max(time)-500])
    ylim(z(1),[SDFyLimMin SDFyLimMax])
    set(sdfFig, 'Position', get(0,'Screensize'))
    str = [sessionName ' ' cellName];
    str = regexprep(str, '\_', ' ');
    suptitle(str)
    
    
    if p.Results.saveFigFlag == 1 
        saveFigurePDF(sdfFig,[savePath sessionName '_' spikeFields{cellInd(i)} '_Qsdfs'])
    end
    
end