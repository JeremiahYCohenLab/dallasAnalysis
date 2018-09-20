function mdlStruct = rwdHistMaxLHspikes_opMD(sessionName, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('intanFlag', 0);
p.addParameter('revForFlag', 0);
p.addParameter('cellName', ['all']);
p.addParameter('trialList', []);
p.addParameter('modelList', {'oneAlpha', 'twoAlpha', 'pe'});
p.addParameter('startValueCSVs', {'spikeModel_startValues.csv', 'spikeModel_2alpha_startValues.csv', 'spikeModel_pe_startValues.csv'});
p.addParameter('spikeFeat', 'preCSspikeCount');
p.parse(varargin{:});

cellName = p.Results.cellName;

%set ranges for model params
alpha_range = [0 1];
tau_range = [0 256];
slope_range = [1 100];
intercept_range = [0 30];

% set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

% define file path
[root,sep] = currComputer();

[animalName] = strtok(sessionName, 'd');
animalName = animalName(2:end);

if isstrprop(sessionName(end), 'alpha')
    sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
else
    sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
end
sortedFolder = dir(sortedFolderLocation);

% get behavior and ephys data
if p.Results.intanFlag
    if any(~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'))) == 1
        sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'));
        load([sortedFolderLocation sortedFolder(sessionDataInd).name])
    else
        [sessionData] = generateSessionData_intan_operantMatching(sessionName);
    end
    [s] = behAnalysisNoPlot_opMD(sessionName, 'revForFlag', p.Results.revForFlag);
    rewardIndsL = find([s.behSessionData.rewardL] == 1); rewardIndsR = find([s.behSessionData.rewardR] == 1);
    for i = 1:length(rewardIndsL)
        sessionData(rewardIndsL(i)).rewardL = 1;
    end
    for i = 1:length(rewardIndsR)
        s.behSessionData(rewardIndsR(i)).rewardR = 1;
    end
else
    if any(~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'))) == 1
        sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'));
        load([sortedFolderLocation sortedFolder(sessionDataInd).name])
    else
        [sessionData] = generateSessionData_nL_operantMatching(sessionName);
    end
end

Qmdls = qLearning_fitOpponency([sessionName '.asc'], p.Results.revForFlag);
pe = Qmdls.fourParams_twoLearnRates_alphaForget.pe;


%% get relevant behavioral data

responseInds = find(~isnan([sessionData.rewardTime]));
outcome = zeros(1, length(responseInds));
allReward_L = [sessionData(responseInds).rewardL];
allReward_R = [sessionData(responseInds).rewardR]; 
outcome(find(allReward_L == 1)) = 1;                
outcome(find(allReward_R == 1)) = 1;    % generate binary outcomes array by trial

choiceTimes = ceil(([sessionData(responseInds).rewardTime] - sessionData(1).CSon)/1000); %baseline to start time and convert to s from ms
if choiceTimes(1) == 0
    choiceTimes(1) = 1;
end
rwdTimes = choiceTimes(logical(outcome));

sessionTime = ceil((sessionData(end).CSon + 3000 - sessionData(1).CSon)/1000);     % find total session time and pad time for reward on last trial
sessionRwds = zeros(1,sessionTime);
sessionRwds(rwdTimes) = 1;      %binary outcomes array by time

iti = diff([sessionData(responseInds).CSon]) / 1000;  %get iti lengths in terms of s


%% extract spike data and put into matrix

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
end

% sometimes no licks/spikes are considered 1x0 and sometimes they are []
% plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
allTrial_spike(cellfun(@isempty,allTrial_spike)) = {zeros(1,0)}; 


%% set time window and smoothing parameters, run analysis for all cells

smoothWin = 250;
trialBeg = tb*1000;
CSoff = tb*1000 + 500;
for i = 1:length(sessionData)
    trialDurDiff(i) = (sessionData(i).trialEnd - sessionData(i).CSon)- tf*1000;
end

mdlStruct = struct;


%% main loop for analyzing each cell

for i = 1:length(cellInd)
    
    if length(cellInd) > 1
        cellNameTemp = cellName{i};
    else
        cellNameTemp = cellName;
    end
    
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
    
    % get relevant spike counts and rates on each trial
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
    
      
    if regexp(p.Results.spikeFeat, 'preCSspikeCount')
        preCSspikeCount = [0 preCSspikeCount(responseInds(2:end))];
        postCSspikeCount = [0 postCSspikeCount(responseInds(2:end))];
        maxFRcs = [0 round(maxFRcs(responseInds(2:end)))];
        minFRcs = [0 minFRcs(responseInds(2:end))];
        maxFRtime = [0 maxFRtime(responseInds(2:end))];
    else
        preCSspikeCount = preCSspikeCount(responseInds);
        postCSspikeCount = postCSspikeCount(responseInds);
        maxFRcs = round(maxFRcs(responseInds));
        minFRcs = minFRcs(responseInds);
        maxFRtime = maxFRtime(responseInds);
    end
    spikeFeat = eval(p.Results.spikeFeat);  
    % run models
    
    for currMod = 1:length(p.Results.modelList)
        
        % set up parameters and matrices for model
        startValues = csvread(p.Results.startValueCSVs{currMod});
        runs = size(startValues, 1);
        allParams = zeros(size(startValues, 1), size(startValues, 2));
        LH = zeros(size(startValues, 1), 1);
        exitFl = zeros(size(startValues, 1), 1);
        numParam = size(startValues, 2);
        A = [eye(size(startValues, 2)); -eye(size(startValues, 2))];

    
        if strcmp(p.Results.modelList{currMod}, 'oneAlpha')
            b = [alpha_range(2); tau_range(2); slope_range(2); intercept_range(2);...
                -alpha_range(1); -tau_range(1); -slope_range(1); -intercept_range(1)];
            parfor r = 1:runs
                [allParams(r,:), LH(r, :), exitFl(r, :)] = ...
                    fmincon(@spikeModel, startValues(r,:), A, b, [], [], [], [], [], options, outcome, iti, spikeFeat);
            end
            [~,bestFit] = min(LH);
            mdlStruct.(cellNameTemp).oneAlpha.bestParams = allParams(bestFit, :);
            [mdlStruct.(cellNameTemp).oneAlpha.LH, mdlStruct.(cellNameTemp).oneAlpha.frProb, mdlStruct.(cellNameTemp).oneAlpha.frEstimate,...
                mdlStruct.(cellNameTemp).oneAlpha.rwdHx] = ...
                spikeModel(mdlStruct.(cellNameTemp).oneAlpha.bestParams, outcome, iti, spikeFeat);
            mdlStruct.(cellNameTemp).oneAlpha.LH = mdlStruct.(cellNameTemp).oneAlpha.LH * -1;
        elseif strcmp(p.Results.modelList{currMod}, 'twoAlpha')
             b = [alpha_range(2); alpha_range(2); tau_range(2); slope_range(2); intercept_range(2);...
                -alpha_range(1); -alpha_range(1); -tau_range(1); -slope_range(1); -intercept_range(1)];
            parfor r = 1:runs
                [allParams(r,:), LH(r, :), exitFl(r, :)] = ...
                    fmincon(@spikeModel_twoAlpha, startValues(r,:), A, b, [], [], [], [], [], options, outcome, iti, spikeFeat);
            end
            [~,bestFit] = min(LH);
            mdlStruct.(cellNameTemp).twoAlpha.bestParams = allParams(bestFit, :);
            [mdlStruct.(cellNameTemp).twoAlpha.LH, mdlStruct.(cellNameTemp).twoAlpha.frProb, mdlStruct.(cellNameTemp).twoAlpha.frEstimate,...
                mdlStruct.(cellNameTemp).twoAlpha.rwdHx] = ...
                spikeModel_twoAlpha(mdlStruct.(cellNameTemp).twoAlpha.bestParams, outcome, iti, spikeFeat);
            mdlStruct.(cellNameTemp).twoAlpha.LH = mdlStruct.(cellNameTemp).twoAlpha.LH * -1;
        elseif strcmp(p.Results.modelList{currMod}, 'pe')
            b = [alpha_range(2); slope_range(2); intercept_range(2);...
                -alpha_range(1); -slope_range(1); -intercept_range(1)];
            parfor r = 1:runs
                [allParams(r,:), LH(r, :), exitFl(r, :)] = ...
                    fmincon(@spikeModel_pe, startValues(r,:), A, b, [], [], [], [], [], options, outcome, iti, spikeFeat, pe);
            end
            [~,bestFit] = min(LH);
            mdlStruct.(cellNameTemp).pe.bestParams = allParams(bestFit, :);
            [mdlStruct.(cellNameTemp).pe.LH, mdlStruct.(cellNameTemp).pe.frProb, mdlStruct.(cellNameTemp).pe.frEstimate,...
                mdlStruct.(cellNameTemp).pe.peHx] = ...
                spikeModel_pe(mdlStruct.(cellNameTemp).pe.bestParams, outcome, iti, spikeFeat, pe);
            mdlStruct.(cellNameTemp).pe.LH = mdlStruct.(cellNameTemp).pe.LH * -1;
        end
    end
    
   
end

mdlStruct.outcome = outcome;
mdlStruct.iti = iti;
mdlStruct.maxFRcs = maxFRcs;
mdlStruct.postCSspikeCount = postCSspikeCount;
mdlStruct.preCSspikeCount = preCSspikeCount;

figure; hold on;
mag = [0.9 0 1];
teal = [0 1 0.9];
set(gcf,'defaultAxesColorOrder',[teal; mag]);

for i = 1:length(p.Results.modelList)
    subplot(length(p.Results.modelList),1,i)
    yyaxis right; plot(smooth(mdlStruct.(p.Results.spikeFeat),4), 'Color', mag, 'linewidth', 1.4);
    ylabel('spike feature')
    yl = ylim;
    yyaxis left; plot(mdlStruct.(cellNameTemp).(p.Results.modelList{i}).frEstimate, 'Color', teal, 'linewidth', 1.4);
    ylabel('model estimate')
    ylim(yl);
    title(p.Results.modelList{i})
    xlim([1 length(outcome)])
end


cellName = erase(cellNameTemp, '_');
suptitle([sessionName ': ' cellName ' -- ' p.Results.spikeFeat]) 

end



%%
function [LH, frProb, frEstimate, rwdHx] = spikeModel(startValues, outcome, iti, spikeFeat)

alpha = startValues(1);
tau = startValues(2);
slope = startValues(3);
intercept = startValues(4);

rwdHx = zeros(1,length(outcome));

for t = 1:length(outcome) - 1
    rwdHx(t+1) = rwdHx(t) + alpha * (outcome(t));
%    rwdHx(t+1) = rwdHx(t) + alpha * (outcome(t) - rwdHx(t));
    rwdHx(t+1) = rwdHx(t+1) * exp(-iti(t) / tau);
end

frEstimate = rwdHx * slope + intercept;
zeroInds = find(frEstimate < 0);
frProb = poisspdf(spikeFeat, frEstimate);
%frProb = normpdf(spikeFeat, frEstimate);
frProb(zeroInds) = 0;


if any(isinf(log(frProb)))
    LH = 1e9;
else
    LH = -1 * sum(log(frProb));
end

end


%%
function [LH, frProb, frEstimate, rwdHx] = spikeModel_twoAlpha(startValues, outcome, iti, spikeFeat)

alphaR = startValues(1);
alphaNR = startValues(2);
tau = startValues(3);
slope = startValues(4);
intercept = startValues(5);

rwdHx = zeros(1,length(outcome));

for t = 1:length(outcome) - 1
    if outcome(t) == 1
        rwdHx(t+1) = rwdHx(t) + alphaR * (outcome(t));
%        rwdHx(t+1) = rwdHx(t) + alphaR * (outcome(t) - rwdHx(t));
    else
        rwdHx(t+1) = rwdHx(t) + alphaNR * (outcome(t));
%        rwdHx(t+1) = rwdHx(t) + alphaNR * (outcome(t) - rwdHx(t));
    end
    rwdHx(t+1) = rwdHx(t+1) * exp(-iti(t) / tau);
end

frEstimate = rwdHx * slope + intercept;
zeroInds = find(frEstimate < 0);
frProb = poisspdf(spikeFeat, frEstimate);
frProb(zeroInds) = 0;
%frProb = normpdf(spikeFeat, frEstimate);

if any(isinf(log(frProb)))
    LH = 1e9;
else
    LH = -1 * sum(log(frProb));
end

end

function [LH, frProb, frEstimate, peHx] = spikeModel_pe(startValues, outcome, iti, spikeFeat, pe)

alpha = startValues(1);
slope = startValues(2);
intercept = startValues(3);

peHx = zeros(1,length(pe));

for t = 1:length(pe) - 1
    peHx(t+1) = peHx(t) + alpha * (abs(pe(t)) - peHx(t));
end

peHx = peHx / max(peHx);
frEstimate = peHx * slope + intercept;
zeroInds = find(frEstimate < 0);
frProb = poisspdf(spikeFeat, frEstimate);
frProb(zeroInds) = 0;


if any(isinf(log(frProb)))
    LH = 1e9;
else
    LH = -1 * sum(log(frProb));
end

end