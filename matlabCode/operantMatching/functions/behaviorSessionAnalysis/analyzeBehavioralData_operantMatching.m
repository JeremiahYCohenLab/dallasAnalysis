function [allRewards, allChoices, glm_rwd, glm_choice, glm_rwdANDchoice, choiceSlope, rwdSlope, lickLat_short, blockSwitch, blockProbs] = analyzeBehavioralData_operantMatching(fileOrFolder, saveFigFlag)

% [allRewards, allChoices, glm_rwd, glm_choice, glm_rwdANDchoice] = analyzeBehavioralData_operantMatching(fileOrFolder, saveFigFlag)
% INPUTS
%   fileOrFolder: str of either .asc file (with .asc extension) or home folder (e.g. mBB040d20161018)
%   saveFigFlag: 1 by default
% OUTPUTS
%   allRewards: binary vector of rewards
%   allChoices: vector of choices; -1 left, 1 right
%   lm_...: linear regression models


if nargin < 2
    saveFigFlag = 1;
end

% % Add current folder to working directory
% currPath = pwd;
% tmp = strfind(currPath,'\');
% addpath(genpath(currPath(1:tmp(end)))); 

[root, sep] = currComputer();

if ~isempty(strfind(fileOrFolder,'.asc')) % input is .asc file
    filename = fileOrFolder;
    [animalName, date] = strtok(filename, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];
    behavioralDataPath = [root animalName sep sessionFolder sep 'behavior' sep filename];
    sessionName = filename(1:strfind(filename,'.asc')-1);
    saveFigName = sessionName;
else % input is the folder
    sessionFolder = fileOrFolder;
    animalName = strtok(sessionFolder, 'd');
    animalName = animalName(2:end);
    filepath = [root animalName sep sessionFolder sep 'behavior' sep];
    allFiles = dir(filepath);
    fileInd = ~cellfun(@isempty,strfind({allFiles.name},'.asc'));
    behavioralDataPath = [filepath allFiles(fileInd).name];
    sessionName = allFiles(fileInd).name;
    saveFigName = sessionFolder(~(sessionFolder==sep));
end

sortedFolderLocation = [root animalName sep sessionFolder sep 'sorted' sep];
sortedFolder = dir(sortedFolderLocation);
if any(~cellfun(@isempty,strfind({sortedFolder.name},[sessionName '_behav.mat'])))
    sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_behav.mat'));
    load([sortedFolderLocation sortedFolder(sessionDataInd).name])
else
    [sessionData, blockSwitch, blockProbs] = generateSessionData_behav_operantMatching(sessionName);
end


%% Break session down into CS+ trials where animal responded

responseInds = find(~isnan([sessionData.rewardTime]));
omitInds = isnan([sessionData.rewardTime]);

origBlockSwitch = blockSwitch; % because I change blockSwitch below and need original below
tempBlockSwitch = blockSwitch;
for i = 2:length(blockSwitch)
    subVal = sum(omitInds(tempBlockSwitch(i-1):tempBlockSwitch(i)));
    blockSwitch(i:end) = blockSwitch(i:end) - subVal;
end

allReward_R = [sessionData(responseInds).rewardR]; 
allReward_L = [sessionData(responseInds).rewardL]; 
allChoices = NaN(1,length(sessionData(responseInds)));
allChoices(~isnan(allReward_R)) = 1;
allChoices(~isnan(allReward_L)) = -1;

allReward_R(isnan(allReward_R)) = 0;
allReward_L(isnan(allReward_L)) = 0;
allChoice_R = double(allChoices == 1);
allChoice_L = double(allChoices == -1);

allRewards = zeros(1,length(allChoices));
allRewards(logical(allReward_R)) = 1;
allRewards(logical(allReward_L)) = -1;

%% 
figure
set(gcf, 'Position', get(0,'Screensize'))
suptitle(sessionName)


%%
subplot(6,8,[17:20 25:28]); hold on
normKern = normpdf(-15:15,0,4);
normKern = normKern / sum(normKern);
xVals = (1:(length(normKern) + length(allChoices) - 1)) - round(length(normKern)/2);
plot(xVals, conv(allChoices,normKern)/max(conv(allChoices,normKern)),'k','linewidth',2);
plot(xVals, conv(allRewards,normKern)/max(conv(allRewards,normKern)),'--','Color',[100 100 100]./255,'linewidth',2)
xlabel('Trials')
ylabel('<-- Left       Right -->')
legend('Choices','Rewards')
xlim([1 length(allChoice_R)])
ylim([-1 1])

for i = 1:length(blockSwitch)
    bs_loc = blockSwitch(i);
    plot([bs_loc bs_loc],[-1 1],'--','linewidth',1,'Color',[30 144 255]./255)
    text(bs_loc,1.04,num2str(blockProbs{i}));
    set(text,'FontSize',3);
end
text(0,1.12,'L/R');

%% linear regression model
tMax = 10;
rwdMatx = [];
for i = 1:tMax
    rwdMatx(i,:) = [NaN(1,i) allRewards(1:end-i)];
end

choiceMatx = [];
for i = 1:tMax
    choiceMatx(i,:) = [NaN(1,i) allChoices(1:end-i)];
end


glm_rwd = fitglm([rwdMatx]', allChoice_R,'distribution','binomial','link','logit'); rsq{1} = num2str(round(glm_rwd.Rsquared.Adjusted*100)/100);
glm_choice = fitglm([choiceMatx]', allChoice_R, 'distribution','binomial','link','logit'); rsq{3} = num2str(round(glm_choice.Rsquared.Adjusted*100)/100);
glm_rwdANDchoice = fitglm([rwdMatx; choiceMatx]', allChoice_R, 'distribution','binomial','link','logit'); rsq{2} = num2str(round(glm_rwdANDchoice.Rsquared.Adjusted*100)/100);

% CIbands_rwd = coefCI(glm_rwd);
% CIbands_rwdANDchoice = coefCI(glm_rwdANDchoice);

subplot(6,8,[39 40]); hold on
relevInds = 2:tMax+1;
coefVals = glm_rwd.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwd);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:tMax,coefVals,errorL,errorU,'k','linewidth',2)

relevInds = 2:tMax+1;
coefVals = glm_rwdANDchoice.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwdANDchoice);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[189,183,107]/255,'linewidth',2)
xlabel('Reward n Trials Back')
ylabel('\beta Coefficient')
% legend(['R-R^2:' rsq{1}],['R&C-R^2:' rsq{2}],'location','northeast')
legend(sprintf('%s\n%s',['R-R^s:' rsq{1}], ['Int: ' num2str(round(100*glm_rwd.Coefficients.Estimate(1))/100)]), ...
       sprintf('%s\n%s\n%s',['R&C-R^2' rsq{2}], ['Int: ' num2str(round(100*glm_rwdANDchoice.Coefficients.Estimate(1))/100)]), ...
       'location','northeast')
xlim([0.5 tMax+0.5])
plot([0 tMax],[0 0],'k--')


subplot(6,8,[47 48]); hold on
relevInds = 2:tMax+1;
coefVals = glm_choice.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_choice);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:tMax,coefVals,errorL,errorU,'Color',[34,139,34]/255,'linewidth',2)

% glm_rwdANDchoice = fitlm([rwdMatx; choiceMatx]', allChoices, 'linear');
relevInds = tMax+2:length(glm_rwdANDchoice.Coefficients.Estimate);
coefVals = glm_rwdANDchoice.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwdANDchoice);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[189,183,107]/255,'linewidth',2)
xlabel('Choice n Trials Back')
ylabel('\beta Coefficient')
legend(['C-R^2:' rsq{3}],['R&C-R^2:' rsq{2}],'location','northeast');
xlim([0.5 tMax+0.5])
plot([0 tMax],[0 0],'k--')

%%
changeChoice = [false abs(diff(allChoices)) > 0];
changeHistogram = [];
changeHistogram_LtoR = [];
changeHistogram_RtoL = [];
changeTimeHistogram_LtoR = [];
changeTimeHistogram_RtoL = [];
for i = find(changeChoice == 1)
    if allChoices(i) == 1 % if a right lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allChoices(i-goBack) == -1 && allRewards(i-goBack) == 0 % if previous trial was a L lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram = [changeHistogram temp];
        changeHistogram_LtoR = [changeHistogram_LtoR temp];
        if i - goBack > 0
            changeTimeHistogram_LtoR = [changeTimeHistogram_LtoR (sessionData(responseInds(i)).CSon - sessionData(responseInds(i-goBack)).rewardTime)];
        end
    elseif allChoices(i) == -1 %if a left lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allChoices(i-goBack) == 1 && allRewards(i-goBack) == 0 % if previous trial was a R lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram = [changeHistogram temp];
        changeHistogram_RtoL = [changeHistogram_RtoL temp];
        if i - goBack > 0
            changeTimeHistogram_RtoL = [changeTimeHistogram_RtoL (sessionData(responseInds(i)).CSon - sessionData(responseInds(i-goBack)).rewardTime)];
        end
    end
end

subplot(6,8,[31 32]); hold on;
histogram(changeHistogram_LtoR,0:max(changeHistogram))
histogram(changeHistogram_RtoL,0:max(changeHistogram))
legend('L -> R','R -> L')
xlabel('Consec No Rewards Before Switch')

changeTimeHistogram_LtoR = ceil(changeTimeHistogram_LtoR / 1000);
changeTimeHistogram_RtoL = ceil(changeTimeHistogram_RtoL / 1000);

subplot(6,8,[29 30]); hold on;
histogram(changeTimeHistogram_LtoR,0:10:(max([max(changeTimeHistogram_LtoR) max(changeTimeHistogram_RtoL)])+10))
histogram(changeTimeHistogram_RtoL,0:10:(max([max(changeTimeHistogram_LtoR) max(changeTimeHistogram_RtoL)])+10))
legend('L -> R','R -> L')
xlabel('Time Since Last Reward Before Switch')

%%
subplot(6,8,[23 24]);
stayDuration = diff([1 find(diff(allChoices) ~= 0)]);
if ~isempty(stayDuration)
    histogram(stayDuration,1:max(stayDuration),'Normalization','probability')
    xl = xlim; yl = ylim;
    allRewardsBin = allRewards;
    allRewardsBin(allRewards==-1) = 1;
    changeChoice = changeChoice(2:end);
    probSwitchNoRwd = sum(changeChoice(allRewardsBin(1:end-1)==0))/sum(allRewardsBin(1:end-1)==0);
    probStayRwd = 1 - (sum(changeChoice(allRewardsBin(1:end-1)==1))/sum(allRewardsBin(1:end-1)==1));
    text('Position', [xl(2) yl(2) 0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'string', ...
        sprintf(['P(stay|rwd) = ' num2str(probStayRwd) '\nP(switch|no rwd) = ' num2str(probSwitchNoRwd)]));
end
xlabel('Stay Duration')


%%
cumsum_allRchoice = cumsum(allChoices == 1);
cumsum_allLchoice = cumsum(allChoices == -1);
cumsum_allRreward = cumsum(allRewards == 1);
cumsum_allLreward = cumsum(allRewards == -1);
cumsum_blockSwitch = cumsum_allLchoice(blockSwitch);
if cumsum_blockSwitch(1) == 0
    cumsum_blockSwitch(1) = 1;
end

halfKern = normKern(round(length(normKern)/2):end);
choiceSlope = atand(diff(conv(cumsum_allRchoice,halfKern))./diff(conv(cumsum_allLchoice,halfKern)));
rwdSlope = atand(diff(conv(cumsum_allRreward,halfKern))./diff(conv(cumsum_allLreward,halfKern)));

subplot(6,8,[37 38 45 46]); hold on
plot(cumsum_allLchoice, cumsum_allRchoice,'linewidth',2,'Color',[30,144,255]/255);

avgRwdSlope = [];
tempMax = 0;
for i = 1:length(cumsum_blockSwitch)
    if i ~= length(cumsum_blockSwitch)
        avgRwdSlope(i) = tand(nanmean(rwdSlope(blockSwitch(i):blockSwitch(i+1))));
        xval = [cumsum_blockSwitch(i) cumsum_blockSwitch(i+1)];
        yval = [cumsum_blockSwitch(i) cumsum_blockSwitch(i+1)]*avgRwdSlope(i) - cumsum_blockSwitch(i)*avgRwdSlope(i) + cumsum_allRchoice(blockSwitch(i));
        tempMax = yval(2);
        plot(xval, yval, 'k','linewidth',2);
    else
        avgRwdSlope(i) = tand(mean(rwdSlope(blockSwitch(i):end)));
        xval = [cumsum_blockSwitch(i) cumsum_allLchoice(end)];
        yval = [cumsum_blockSwitch(i) cumsum_allLchoice(end)]*avgRwdSlope(i) - cumsum_blockSwitch(i)*avgRwdSlope(i) + cumsum_allRchoice(blockSwitch(i));
        tempMax = yval(2);
        plot(xval, yval, 'k','linewidth',2);
    end
end
limMax = max([max(cumsum_allLchoice) max(cumsum_allRchoice)]);
xlim([0 limMax])
ylim([0 limMax])
legend('Choice','Income','location','best')
xlabel('Cumulative Left Choices'); ylabel('Cumulative Right Choices')

%%
subplot(6,8,[33:36 41:44]); hold on;

plot(choiceSlope,'linewidth',1.5,'Color',[30,144,255]/255)
plot(rwdSlope,'k','linewidth',1.5)
for i = 1:length(avgRwdSlope)
    if i ~= length(avgRwdSlope)
        plot([blockSwitch(i) blockSwitch(i+1)], [mean(choiceSlope(blockSwitch(i):blockSwitch(i+1))) mean(choiceSlope(blockSwitch(i):blockSwitch(i+1)))],'linewidth',4,'Color',[30,144,255]/255);
        plot([blockSwitch(i) blockSwitch(i+1)], [mean(rwdSlope(blockSwitch(i):blockSwitch(i+1))) mean(rwdSlope(blockSwitch(i):blockSwitch(i+1)))],'k','linewidth',4);
    else
        plot([blockSwitch(i) length(choiceSlope)], [mean(choiceSlope(blockSwitch(i):length(choiceSlope))) mean(choiceSlope(blockSwitch(i):length(choiceSlope)))],'linewidth',4,'Color',[30,144,255]/255)
        plot([blockSwitch(i) length(rwdSlope)], [mean(rwdSlope(blockSwitch(i):length(choiceSlope))) mean(rwdSlope(blockSwitch(i):length(choiceSlope)))],'k','linewidth',4);
    end
end

xlim([1 length(allChoices)])
ylim([0 90])
legend('Choice','Rewards')
xlabel('Choice number')
ylabel('Slope (degrees)')


%%
rMag = 1;
nrMag = rMag/2;

% trial plot
subplot(6,8,[1:8]); hold on
xlabel('Trials')
xlim([1 length(sessionData)])
ylim([-1 1])
for i = 1:length(sessionData)
    if strcmp(sessionData(i).trialType,'CSplus')
        if ~isnan(sessionData(i).rewardR)
            if sessionData(i).rewardR == 1 % R side rewarded
                plot([i i],[0 rMag],'k')
            else
                plot([i i],[0 nrMag],'k') % R side not rewarded
            end
        elseif ~isnan(sessionData(i).rewardL)
            if sessionData(i).rewardL == 1 % L side rewarded
                plot([i i],[-1*rMag 0],'k')
            else
                plot([i i],[-1*nrMag 0],'k')
            end
            
        else % CSplus trial but no rewardL or rewardR
            plot([i i],[-rMag rMag],'r')
        end
    else % CS minus trial
        plot([i],0,'ko','markersize',4,'linewidth',2)
        if ~isempty(sessionData(i).licksL) & sessionData(i).licksL(1) - sessionData(i).CSon < 2000
            plot([i i],[-1*nrMag 0],'k')
        elseif ~isempty(sessionData(i).licksR) & sessionData(i).licksR(1) - sessionData(i).CSon < 2000
            plot([i i],[nrMag 0],'k')
        end 
    end
    if any(i == origBlockSwitch)
        plot([i i],[-1*rMag rMag],'--','linewidth',1,'Color',[30 144 255]./255)
    end
end

for i = 1:length(blockSwitch)
    bs_loc = origBlockSwitch(i);
    plot([bs_loc bs_loc],[-1 1],'--','linewidth',1,'Color',[30 144 255]./255)
    text(bs_loc,1.12,blockProbs{i});
    set(text,'FontSize',3);
end
text(0,1.35,'L/R');
ylabel('<-- L       R  -->')


% time plot
subplot(6,8,[9:16]); hold on
xlabel('Time (min)')
for i = 1:length(sessionData)
    currTime = (sessionData(i).CSon - sessionData(1).CSon)/1000/60; %convert to min
    if strcmp(sessionData(i).trialType,'CSplus')
        if ~isnan(sessionData(i).rewardR)
            if sessionData(i).rewardR == 1 % R side rewarded
                plot([currTime currTime],[0 rMag],'k')
            else
                plot([currTime currTime],[0 nrMag],'k') % R side not rewarded
            end
        elseif ~isnan(sessionData(i).rewardL)
            if sessionData(i).rewardL == 1 % L side rewarded
                plot([currTime currTime],[-1*rMag 0],'k')
            else
                plot([currTime currTime],[-1*nrMag 0],'k')
            end
            
        else % CSplus trial but no rewardL or rewardR
            plot([currTime currTime],[-rMag rMag],'r')
        end
    else % CS minus trial
        plot([currTime],0,'ko','markersize',4,'linewidth',2)
        if ~isempty(sessionData(i).licksL) & sessionData(i).licksL(1) - sessionData(i).CSon < 2000
            plot([currTime currTime],[-1*nrMag 0],'k')
        elseif ~isempty(sessionData(i).licksR) & sessionData(i).licksR(1) - sessionData(i).CSon < 2000
            plot([currTime currTime],[nrMag 0],'k')
        end 
    end
    if any(i == origBlockSwitch)
        plot([currTime currTime],[-1*rMag rMag],'--','linewidth',1,'Color',[30 144 255]./255)
    end
end
xlim([0 currTime]);

%% determine and plot lick latency distributions for each spout
lickLat = [];       lickRate = [];
lickLat_L = [];     lickRate_L = [];
lickLat_R = [];     lickRate_R = [];
for i = 1:length(sessionData)
    if ~isempty(sessionData(i).rewardTime)
        lickLat = [lickLat sessionData(i).rewardTime - sessionData(i).CSon];
        if ~isnan(sessionData(i).rewardL)
            lickLat_L = [lickLat_L sessionData(i).rewardTime - sessionData(i).CSon];
            if sessionData(i).rewardL == 1
                if length(sessionData(i).licksL) > 1
                    lickRateTemp = 1000/(min(diff(sessionData(i).licksL)));
                    lickRate = [lickRate lickRateTemp];
                    lickRate_L = [lickRate_L lickRateTemp];
                else
                   lickRate = [lickRate 0];
                   lickRate_L = [lickRate_L 0]; 
                end
            end
        elseif ~isnan(sessionData(i).rewardR)
            lickLat_R = [lickLat_R sessionData(i).rewardTime - sessionData(i).CSon];      %make single licks zeros for easier indexing
            if sessionData(i).rewardR == 1
                if length(sessionData(i).licksR) > 1
                    lickRateTemp = 1000/(min(diff(sessionData(i).licksR)));
                    lickRate = [lickRate lickRateTemp];
                    lickRate_R = [lickRate_R lickRateTemp];
                else
                    lickRate = [lickRate 0];
                    lickRate_R = [lickRate_R 0];
                end
            end
        end
    end
end

subplot(6,8,[21 22]); hold on
histogram(lickLat_L,0:50:1500,'Normalization','probability'); histogram(lickLat_R,0:50:1500,'Normalization','probability')
legend('Left Licks','Right Licks')
xlabel('Lick Latency (ms)')

%% Z-scored lick latency analysis

lickLatResp = lickLat(responseInds);                    %remove NaNs from lickLat array
lickLatResp = lickLatResp(2:end);                       %shift for comparison to rwd history
lickLatInds = find(lickLatResp > 250);                  %find indices of non-preemptive licks (limit to normal distribution)

if ~isnan(sessionData(responseInds(1)).rewardR)      %remove first response for shift to compare to rwd hist
    responseLat_R = lickLat_R(2:end);
    responseLat_L = lickLat_L;
else
    responseLat_R = lickLat_R;
    responseLat_L = lickLat_L(2:end);
end

responseLat_R = responseLat_R(responseLat_R > 250);        %remove lick latencies outside of normal distribution
responseLat_L = responseLat_L(responseLat_L > 250);
responseLat_R  = zscore(responseLat_R);                   %get z scores for lick latencies based on spout side average
responseLat_L  = zscore(responseLat_L);
choicesLick = allChoices(2:end);                        %make shifted choice array without preemptive licks
choicesLick = choicesLick(lickLatInds);

L = 1;
R = 1;
for j = 1:length(choicesLick)                     %put z scored lick latencies back in trial order
    if choicesLick(j) == 1
        responseLat(j) = responseLat_R(R);
        R = R + 1;
    else
        responseLat(j) = responseLat_L(L);
        L = L + 1;
    end
end

respRange = minmax(responseLat);
if -respRange(1,1) > respRange(1,2)
    normRespLat = responseLat / -respRange(1,1);
else
    normRespLat = responseLat / respRange(1,2);
end

%%

if saveFigFlag == 1
    savepath = [root animalName sep sessionFolder sep  'figures' sep];
    if isempty(dir(savepath))
        mkdir(savepath)
    end

    % sessionFolder(~(sessionFolder==sep)) removes the separator (/ or \) in the filename before saving
    saveFigurePDF(gcf,[savepath sep saveFigName '_behavior'])
end

%%
figure   %make new lick behavior analysis figure
suptitle(sessionName)


%% lick latency and recent rwd hist analysis

timeMax = 61000;
binSize = 10000;
timeBinEdges = [1000:binSize:timeMax];  %no trials shorter than 1s between outcome and CS on
tMax = length(timeBinEdges) - 1;
rwdTimeMatx = zeros(tMax, length(responseInds));     %initialize matrices for number of response trials x number of time bins
for j = 2:length(responseInds)          
    k = 1;
    %find time between "current" choice and previous rewards, up to timeMax in the past 
    timeTmp = []; timeTmpN =[];
    while j-k > 0 & sessionData(responseInds(j)).rewardTime - sessionData(responseInds(j-k)).rewardTime < timeMax
        if sessionData(responseInds(j-k)).rewardL == 1 || sessionData(responseInds(j-k)).rewardR == 1
            timeTmp = [timeTmp (sessionData(responseInds(j)).rewardTime - sessionData(responseInds(j-k)).rewardTime)];
        end
        k = k + 1;
    end
    %bin outcome times and use to fill matrices
    if ~isempty(timeTmp)
        binnedRwds = discretize(timeTmp,timeBinEdges);
        for k = 1:tMax
            if ~isempty(find(binnedRwds == k))
                rwdTimeMatx(k,j) = sum(binnedRwds == k);
            end
        end
    end
end

%fill in NaNs at beginning of session
j = 2;
while sessionData(responseInds(j)).rewardTime - sessionData(responseInds(1)).rewardTime < timeMax
    tmpDiff = sessionData(responseInds(j)).rewardTime - sessionData(responseInds(1)).rewardTime;
    binnedDiff = discretize(tmpDiff, timeBinEdges);
    rwdTimeMatx(binnedDiff:tMax,j) = NaN;
    j = j+1;
end

rwdTimeMatx(:,1) = NaN;
glm_rwdLickLat = fitlm([rwdTimeMatx(:,lickLatInds)]', responseLat);                 
glm_rwdLickRate = fitlm([rwdTimeMatx(:,logical(allRewardsBin))]', lickRate);

%plot beta coefficients from lrm's
subplot(4,2,1);
relevInds = 2:tMax+1;
coefVals = glm_rwdLickLat.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwdLickLat);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

xlabel('reward n seconds back')
ylabel('\beta Coefficient')
xlim([0 (tMax*binSize/1000 + 5)])
title('LRM: rewards in time on lick latency')

subplot(4,2,2);
coefVals = glm_rwdLickRate.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwdLickRate);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

xlabel('reward n seconds back')
ylabel('\beta Coefficient')
xlim([0 (tMax*binSize/1000 + 5)])
title('LRM: rewards in time on lick rate')

%plot reward history smoothed with kernel derived from beta coefficients (licks) v lick behavior 
purp = [0.5 0 0.8];
blue = [0 0 1];
set(gcf,'defaultAxesColorOrder',[blue; purp]);

rwdsSmooth = fastsmooth(allRewardsBin, 5, 3);

subplot(4,2,3)
yyaxis left; plot(rwdsSmooth(lickLatInds), 'b'); hold on
ylabel('Reward History (from lick LRM)')
yyaxis right; plot(fastsmooth(responseLat,5,3),'-','Color', purp);
ylabel('Response Latency Z-Score')
xlabel('Choice Trials')


subplot(4,2,4)
yyaxis left; plot(rwdsSmooth(logical(allRewardsBin)), 'b'); hold on
ylabel('Reward History (from lick LRM)')
yyaxis right; plot(fastsmooth(lickRate,5,3),'-','Color', purp);
ylabel('Peak Lick Rate Z-Score')
xlabel('Rewarded Trials')

%plot scatters for reward history smoothed with kernel derived from beta coefficients (choices) v lick behavior  
expFit = singleExpFit(glm_rwd.Coefficients.Estimate(2:end));
expConv = expFit.a*exp(-(1/expFit.b)*(1:10));
expConv = expConv./sum(expConv);
rwdHx = conv(allRewardsBin,expConv);
rwdHx = rwdHx(1:end-(length(expConv)-1)); 

rwdHxLick = rwdHx(lickLatInds);                   %remove preemptive lick trials from rwd hist

coeffs = polyfit(rwdHxLick, responseLat, 1);
fittedX = linspace(0, 1, 200);
fittedY = polyval(coeffs, fittedX);
mdl = fitlm(rwdHxLick, responseLat);
rSqr = mdl.Rsquared(1).Ordinary(1);

purp = [0.5 0 0.8];
blue = [0 0 1];
set(gcf,'defaultAxesColorOrder',[blue; purp]);

subplot(4,2,5)
scatter(rwdHxLick, responseLat,'b'); hold on
ylabel('Lick Latency Z Score')
xlabel('Reward History (from choice LRM)')
plot(fittedX, fittedY, '-','Color', purp, 'LineWidth', 3);
legend('Licks', sprintf(['R^2: ' num2str(rSqr)]))

%
rewardsLick = allRewards(allRewards == 1 | allRewards == -1);
if sessionData(responseInds(1)).rewardR == 1      %remove first response for shift to compare to rwd hist
    responseRate_R = lickRate_R(2:end);
    responseRate_L = lickRate_L;
    responseRate = lickRate(2:end);
    responseRateInds = find(lickRate(2:end) < 15);
    rewardsLick = rewardsLick(2:end);
elseif sessionData(responseInds(1)).rewardL == 1 
    responseRate_R = lickRate_R;
    responseRate_L = lickRate_L(2:end);
    responseRate = lickRate(2:end);
    responseRateInds = find(lickRate(2:end) < 15);
    rewardsLick = rewardsLick(2:end);
else
    responseRate_R = lickRate_R;
    responseRate_L = lickRate_L;
    responseRate = lickRate;
    responseRateInds = find(lickRate < 15);
end


responseRate = responseRate(responseRate < 15);
responseRate_R = responseRate_R(responseRate_R < 15);
responseRate_L = responseRate_L(responseRate_L < 15);
responseRate_R = zscore(responseRate_R);
responseRate_L = zscore(responseRate_L);
rewardsLick = rewardsLick(responseRateInds);

L = 1;
R = 1;
for j = 1:length(rewardsLick)                     %put z scored lick rates back in trial order
    if rewardsLick(j) == 1
        responseRate(j) = responseRate_R(R);
        R = R + 1;
    elseif rewardsLick(j) == -1
        responseRate(j) = responseRate_L(L);
        L = L + 1;
    end
end


binAllRewards = logical(allRewards == 1 | allRewards == -1);
rwdHxRwd = rwdHx(binAllRewards(2:end)==1);
rwdHxRwd = rwdHxRwd(responseRateInds);


coeffs2 = polyfit(rwdHxRwd, responseRate, 1);
fittedY2 = polyval(coeffs2, fittedX);
mdl2 = fitlm(rwdHxRwd, responseRate);
rSqr2 = mdl2.Rsquared(1).Ordinary(1);

subplot(4,2,6)
scatter(rwdHxRwd, responseRate, 'b'); hold on;
xlabel('Reward History (from choice LRM)')
ylabel('Peak Lick Rate Z-Score')
plot(fittedX, fittedY2, '-','Color', purp, 'LineWidth', 3);
legend('Licks', sprintf(['R^2: ' num2str(rSqr2)]))


subplot(4,2,7)
yyaxis left; plot(fastsmooth(rwdHxLick, 5, 3), 'b'); hold on
ylabel('Reward History (from choice LRM)')
yyaxis right; plot(fastsmooth(responseLat,5,3),'-','Color', purp);
ylabel('Response Latency Z-Score')
xlabel('Choice Trials')

subplot(4,2,8)
yyaxis left; plot(fastsmooth(rwdHxRwd, 5, 3), 'b'); hold on
ylabel('Reward History (from choice LRM)')
yyaxis right; plot(fastsmooth(responseRate,5,3),'-','Color', purp);
ylabel('Peak Lick Rate Z-Score')
xlabel('Rewarded Trials')


set(gcf, 'Position', get(0,'Screensize'))

savePath = [root animalName sep sessionFolder sep  'figures' sep];
if isempty(dir(savePath))
    mkdir(savePath)
end

if saveFigFlag == 1
    saveFigurePDF(gcf,[savePath sep sessionName '_lickBehavior'])
end
