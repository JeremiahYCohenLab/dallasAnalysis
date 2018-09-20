function [glm_rwd, glm_rwdTime] = compareLogRegResiduals_opMD(sessionName)

[root, sep] = currComputer();

[animalName, date] = strtok(sessionName, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];

if isstrprop(sessionName(end), 'alpha')
    sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
else
    sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
end

if exist(sessionDataPath,'file')
    load(sessionDataPath)
else
    [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName);
end


%create arrays for choices and rewards
responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
omitInds = isnan([behSessionData.rewardTime]); 

origBlockSwitch = blockSwitch;
tempBlockSwitch = blockSwitch;
for i = 2:length(blockSwitch)
    subVal = sum(omitInds(tempBlockSwitch(i-1):tempBlockSwitch(i)));
    blockSwitch(i:end) = blockSwitch(i:end) - subVal;
end

allReward_R = [behSessionData(responseInds).rewardR]; 
allReward_L = [behSessionData(responseInds).rewardL]; 
allChoices = NaN(1,length(behSessionData(responseInds)));
allChoices(~isnan(allReward_R)) = 1;
allChoices(~isnan(allReward_L)) = -1;
allChoice_R = double(allChoices == 1);
allChoice_L = double(allChoices == -1);
changeChoice = [false abs(diff(allChoices)) > 0];

allReward_R(isnan(allReward_R)) = 0;
allReward_L(isnan(allReward_L)) = 0;
allRewards = zeros(1,length(allChoices));
allRewards(logical(allReward_R)) = 1;
allRewards(logical(allReward_L)) = -1;

allProbsL = [behSessionData(responseInds).rewardProbL];
allProbsR = [behSessionData(responseInds).rewardProbR];

tMax = 20;
timeMax = 127000;
binSize = (timeMax - 1000)/tMax;
timeBinEdges = [1000:binSize:timeMax];

rwdMatx = [];
for j = 1:tMax
    rwdMatx(j,:) = [NaN(1,j) allRewards(1:end-j)];
end
glm_rwd = fitglm([rwdMatx]', allChoice_R,'distribution','binomial','link','logit');
aicTrial = aicbic(glm_rwd.LogLikelihood, tMax);


%logistic regression for rwds in time
rwdTimeMatx = zeros(tMax, length(responseInds));     %initialize matrices for number of response trials x number of time bins
noRwdTimeMatx = zeros(tMax, length(responseInds));
for j = 2:length(responseInds)          
    k = 1;
    %find time between "current" choice and previous rewards, up to timeMax in the past 
    timeTmpL = []; timeTmpR = [];
    while j-k > 0 & behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime < timeMax
        if behSessionData(responseInds(j-k)).rewardL == 1
            timeTmpL = [timeTmpL (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
        end
        if behSessionData(responseInds(j-k)).rewardR == 1
            timeTmpR = [timeTmpR (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
        end
        k = k + 1;
    end
    %bin outcome times and use to fill matrices
    if ~isempty(timeTmpL)
        binnedRwds = discretize(timeTmpL,timeBinEdges);
        for k = 1:tMax
            rwdTimeMatx(k,j) = -1*sum(binnedRwds == k);
        end
    end
    if ~isempty(timeTmpR)
        binnedRwds = discretize(timeTmpR,timeBinEdges);
        for k = 1:tMax
            rwdTimeMatx(k,j) = sum(binnedRwds == k);
        end
    end
end

rwdTimeMatx(:,1) = NaN;
glm_rwdTime = fitglm([rwdTimeMatx]', allChoice_R,'distribution','binomial','link','logit'); 
aicTime = aicbic(glm_rwdTime.LogLikelihood, tMax);


%%
figure;
suptitle(sessionName)

%% Plot Raw Data

rMag = 1;
nrMag = rMag/2;

% trial plot
subplot(5,8,[1:8]); hold on
for i = 1:length(responseInds)
    if ~isnan(behSessionData(responseInds(i)).rewardR)
        if behSessionData(responseInds(i)).rewardR == 1 % R side rewarded
            plot([i i],[0 rMag],'k')
        else
            plot([i i],[0 nrMag],'k') % R side not rewarded
        end
    elseif ~isnan(behSessionData(responseInds(i)).rewardL)
        if behSessionData(responseInds(i)).rewardL == 1 % L side rewarded
            plot([i i],[-1*rMag 0],'k')
        else
            plot([i i],[-1*nrMag 0],'k')
        end
    end
end

for i = 1:length(blockSwitch)
    bs_loc = blockSwitch(i);
    if i > 1
        plot([bs_loc bs_loc],[-1 1],'--','linewidth',1,'Color',[30 144 255]./255)
    end
    labelOffset = 1.12;
    a = num2str(allProbsL(blockSwitch(i)+1));
    b = '/';
    c = num2str(allProbsR(blockSwitch(i)+1));
    label = strcat(a,b,c);
    text(bs_loc,labelOffset,label);
    set(text,'FontSize',3);
end
xlabel('Trial')
text(0,1.35,'L/R');
ylabel('<-- L       R  -->')

% time plot
subplot(5,8,[17:24]); hold on
for i = 1:length(responseInds)
    currTime = (behSessionData(responseInds(i)).CSon - behSessionData(1).CSon)/1000/60; %convert to min
    if ~isnan(behSessionData(responseInds(i)).rewardR)
        if behSessionData(responseInds(i)).rewardR == 1 % R side rewarded
            plot([currTime currTime],[0 rMag],'k')
        else
            plot([currTime currTime],[0 nrMag],'k') % R side not rewarded
        end
    elseif ~isnan(behSessionData(responseInds(i)).rewardL)
        if behSessionData(responseInds(i)).rewardL == 1 % L side rewarded
            plot([currTime currTime],[-1*rMag 0],'k')
        else
            plot([currTime currTime],[-1*nrMag 0],'k')
        end
    end
    if any(i == blockSwitch & i~= 1)
        plot([currTime currTime],[-1*rMag rMag],'--','linewidth',1,'Color',[30 144 255]./255)
    end
end
xlim([0 currTime]);
xlabel('Time (min)')

%% plot probability overlayed on trial

subplot(5,8,[9:16]); hold on;
for i = 1:length(responseInds)
    if ~isnan(behSessionData(responseInds(i)).rewardR)
        if behSessionData(i).rewardR == 1 % R side rewarded
            plot([i i],[0.5 1],'k')
        else
            plot([i i],[0.5 0.75],'k') % R side not rewarded
        end
    elseif ~isnan(behSessionData(responseInds(i)).rewardL)
        if behSessionData(responseInds(i)).rewardL == 1 % L side rewarded
            plot([i i],[0 0.5],'k')
        else
            plot([i i],[0.25 0.5],'k')
        end
    end
end

for i = 1:length(blockSwitch)
    if i > 1
        plot([blockSwitch(i) blockSwitch(i)],[-1 1],'--','linewidth',1,'Color',[30 144 255]./255)
    end
end

trialPlot = plot(glm_rwd.Fitted.Probability, 'b','linewidth',2);
timePlot = plot(glm_rwdTime.Fitted.Probability, 'r','linewidth',2);

xlabel('Trial')
ylabel('probability R choice')
legend([trialPlot timePlot],'trial', 'time')
ylim([0 1]);


%% plot probability overlayed on time

subplot(5,8,[25:32]); hold on
CSonTimes = ([behSessionData(responseInds).CSon] - behSessionData(1).CSon)/1000/60;
for i = 1:length(responseInds)
    currTime = CSonTimes(i);
    if ~isnan(behSessionData(responseInds(i)).rewardR)
        if behSessionData(responseInds(i)).rewardR == 1 % R side rewarded
            plot([currTime currTime],[0.5 1],'k')
        else
            plot([currTime currTime],[0.5 0.75],'k') % R side not rewarded
        end
    elseif ~isnan(behSessionData(responseInds(i)).rewardL)
        if behSessionData(responseInds(i)).rewardL == 1 % L side rewarded
            plot([currTime currTime],[0 0.5],'k')
        else
            plot([currTime currTime],[0.25 0.5],'k')
        end
    end
    if any(i == blockSwitch & i~= 1)
        plot([currTime currTime],[-1*rMag rMag],'--','linewidth',1,'Color',[30 144 255]./255)
    end
end
xlim([0 currTime]);
xlabel('Time (min)')


trialPlotTime = plot(CSonTimes, glm_rwd.Fitted.Probability, 'b','linewidth',2);
timePlotTime = plot(CSonTimes, glm_rwdTime.Fitted.Probability, 'r','linewidth',2);

xlabel('Time(min)')
ylabel('probability R choice')
legend([trialPlotTime timePlotTime],'trial', 'time')
ylim([0 1]);
