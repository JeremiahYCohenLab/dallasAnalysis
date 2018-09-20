function [s] = behAnalysisNoPlot_opMD(sessionName, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag',0)
p.addParameter('makeFigFlag', 0)
p.addParameter('tMax', 20)
p.addParameter('timeMax', 121000)
p.addParameter('timeBins', 12)
p.parse(varargin{:});


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

if p.Results.revForFlag
    if exist(sessionDataPath,'file')
        load(sessionDataPath)
        behSessionData = sessionData;
    else
        [behSessionData, blockSwitch, blockProbs] = generateSessionData_behav_operantMatching(sessionName);
    end
else
    if exist(sessionDataPath,'file')
        load(sessionDataPath)
    else
        [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
end

%% Break session down into CS+ trials where animal responded

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
changeChoice = [false abs(diff(allChoices)) > 0];

allReward_R(isnan(allReward_R)) = 0;
allReward_L(isnan(allReward_L)) = 0;
allChoice_R = double(allChoices == 1);
allChoice_L = double(allChoices == -1);

allRewards = zeros(1,length(allChoices));
allRewards(logical(allReward_R)) = 1;
allRewards(logical(allReward_L)) = -1;
allRewardsBinary = allRewards;                      %make all rewards have the same value
allRewardsBinary(find(allRewards==-1)) = 1;
rewardsList =  allRewards(find(allRewards~=0));

if blockSwitch(end) == length(allChoices)
    blockSwitch = blockSwitch(1:end-1);
end


%% determine and plot lick latency distributions for each spout
lickLat = [];       lickRate = [];
lickLat_L = [];     lickRate_L = [];
lickLat_R = [];     lickRate_R = [];
lickSide = NaN(1,length(behSessionData));
for i = 1:length(behSessionData)
    if ~isempty(behSessionData(i).rewardTime)
        lickLat = [lickLat behSessionData(i).rewardTime - behSessionData(i).CSon];
        if ~isnan(behSessionData(i).rewardL)
            lickSide(i) = -1;
            lickLat_L = [lickLat_L behSessionData(i).rewardTime - behSessionData(i).CSon];
            if behSessionData(i).rewardL == 1
                if length(behSessionData(i).licksL) > 1
                    lickRateTemp = 1000/(min(diff(behSessionData(i).licksL)));
                    lickRate = [lickRate lickRateTemp];
                    lickRate_L = [lickRate_L lickRateTemp];
                else
                   lickRate = [lickRate 0];
                   lickRate_L = [lickRate_L 0]; 
                end
            end
        elseif ~isnan(behSessionData(i).rewardR)
            lickSide(i) = 1;
            lickLat_R = [lickLat_R behSessionData(i).rewardTime - behSessionData(i).CSon];      
            if behSessionData(i).rewardR == 1
                if length(behSessionData(i).licksR) > 1     
                    lickRateTemp = 1000/(min(diff(behSessionData(i).licksR)));
                    lickRate = [lickRate lickRateTemp];
                    lickRate_R = [lickRate_R lickRateTemp];  
                else                                                                    %make single licks zeros for easier indexing
                    lickRate = [lickRate 0];
                    lickRate_R = [lickRate_R 0];
                end
            end
        end
    end
end

%% Z-scored lick latency (gets rid of preemptive licks)

lickLatResp = lickLat(responseInds);                    %remove NaNs from lickLat array
lickLatResp = lickLatResp(2:end);                       %shift for comparison to rwd history
lickLatInds = find(lickLatResp > 250);                  %find indices of non-preemptive licks (limit to normal distribution)

if ~isnan(behSessionData(responseInds(1)).rewardR)      %remove first response for shift to compare to rwd hist
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

%% Z scored lick rate (eliminates trials with impossible lick rates)

rewardsLick = allRewards(allRewards == 1 | allRewards == -1);
responseRateInds = find(lickRate < 15);

corrLickRate = lickRate(lickRate < 15);
corrLickRate_R = lickRate_R(lickRate_R < 15);
corrLickRate_L = lickRate_L(lickRate_L < 15);
corrLickRate_R = zscore(corrLickRate_R);
corrLickRate_L = zscore(corrLickRate_L);
rewardsLick = rewardsLick(responseRateInds);

L = 1;
R = 1;
for j = 1:length(rewardsLick)                     %put z scored lick rates back in trial order
    if rewardsLick(j) == 1
        corrLickRate(j) = corrLickRate_R(R);
        R = R + 1;
    elseif rewardsLick(j) == -1
        corrLickRate(j) = corrLickRate_L(L);
        L = L + 1;
    end
end


%% linear regression model by trial

rwdMatx = [];
for i = 1:p.Results.tMax
    rwdMatx(i,:) = [NaN(1,i) allRewards(1:end-i)];
end

choiceMatx = [];
for i = 1:p.Results.tMax
    choiceMatx(i,:) = [NaN(1,i) allChoices(1:end-i)];
end

glm_rwd = fitglm([rwdMatx]', allChoice_R,'distribution','binomial','link','logit'); 
glm_choice = fitglm([choiceMatx]', allChoice_R, 'distribution','binomial','link','logit');
glm_rwdANDchoice = fitglm([rwdMatx; choiceMatx]', allChoice_R, 'distribution','binomial','link','logit');


%% use exponential decay fit to beta values from trial-wise logistic regression to smooth reward hist

%from rwd model
expFit = singleExpFit(glm_rwd.Coefficients.Estimate(2:end));
expConv = expFit.a*exp(-(1/expFit.b)*(1:10));
expConv = expConv./sum(expConv);

rwdHx = conv(allRewardsBinary,expConv);              %convolve with exponential decay to give weighted moving average
rwdHx = rwdHx(1:end-(length(expConv)-1));            %to account for convolution padding
rwdHx = [0 rwdHx(1:end-1)];                          %to make rwdHx for each trial prior to its outcome
rwdHx_L = conv(allReward_L,expConv);                 %same convolution but only with L rewards over all trials
rwdHx_L = rwdHx_L(1:end-(length(expConv)-1));  
rwdHx_R = conv(allReward_R,expConv);
rwdHx_R = rwdHx_R(1:end-(length(expConv)-1));

%from choice model
expFitChoice = singleExpFit(glm_choice.Coefficients.Estimate(2:end));
expConvChoice = expFitChoice.a*exp(-(1/expFitChoice.b)*(1:10));
expConvChoice = expConvChoice./sum(expConvChoice);
allChoices_L = +(allChoices == -1);                                       % + converts to double from logical
allChoices_R = +(allChoices == 1);
choiceHx = conv(allChoices,expConvChoice,'full');
choiceHx = choiceHx(1:end-(length(expConvChoice)-1));
choiceHx = [0 choiceHx(1:end-1)];  


%% linear regression model by time

timeBinSize = (p.Results.timeMax - 1000)/p.Results.timeBins;
timeBinEdges = [1000:timeBinSize:p.Results.timeMax];

%logistic regression for rwds in time
rwdTimeMatx = zeros(p.Results.timeBins, length(responseInds));     %initialize matrices for number of response trials x number of time bins
noRwdTimeMatx =  zeros(p.Results.timeBins, length(responseInds));
for j = 2:length(responseInds)          
    k = 1;
    %find time between "current" choice and previous rewards, up to p.Results.timeMax in the past 
    timeTmpL = []; timeTmpR = []; nTimeTmpL = []; nTimeTmpR = [];
    while j-k > 0 & behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime < p.Results.timeMax
        if behSessionData(responseInds(j-k)).rewardL == 1
            timeTmpL = [timeTmpL (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
        end
        if behSessionData(responseInds(j-k)).rewardR == 1
            timeTmpR = [timeTmpR (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
        end
        if behSessionData(responseInds(j-k)).rewardL == 0
            nTimeTmpL = [nTimeTmpL (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
        end
        if behSessionData(responseInds(j-k)).rewardR == 0
            nTtimeTmpR = [nTimeTmpR (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
        end
        k = k + 1;
    end
    %bin outcome times and use to fill matrices
    if ~isempty(timeTmpL)
        binnedRwds = discretize(timeTmpL,timeBinEdges);
        for k = 1:p.Results.timeBins
            if any(binnedRwds == k)
                rwdTimeMatx(k,j) = -1*sum(binnedRwds == k);
            end
        end
    end
    if ~isempty(timeTmpR)
        binnedRwds = discretize(timeTmpR,timeBinEdges);
        for k = 1:p.Results.timeBins
            if any(binnedRwds == k) & isnan(rwdTimeMatx(k,j))
                rwdTimeMatx(k,j) = sum(binnedRwds == k);
            elseif any(binnedRwds == k) & ~isnan(rwdTimeMatx(k,j))
                rwdTimeMatx(k,j) = rwdTimeMatx(k,j) + sum(binnedRwds == k);
            end
        end
    end
    if ~isempty(nTimeTmpL)
        binnedNoRwds = discretize(nTimeTmpL,timeBinEdges);
        for k = 1:p.Results.timeBins
            if any(binnedNoRwds == k)
                noRwdTimeMatx(k,j) = -1*sum(binnedNoRwds == k);
            end
        end
    end
    if ~isempty(nTimeTmpR)
        binnedNoRwds = discretize(nTimeTmpR,timeBinEdges);
        for k = 1:p.Results.timeBins
            if any(binnedNoRwds == k) & isnan(noRwdTimeMatx(k,j))
                noRwdTimeMatx(k,j) = sum(binnedNoRwds == k);
            elseif any(binnedNoRwds == k) & ~isnan(noRwdTimeMatx(k,j))
                noRwdTimeMatx(k,j) = noRwdTimeMatx(k,j) + sum(binnedNoRwds == k);
            end
        end
    end
end

rwdTimeMatx(:,1) = NaN;
noRwdTimeMatx(:,1) = NaN;
rwdTimeMatxBin = rwdTimeMatx;
rwdTimeMatxBin(find(rwdTimeMatxBin < 0)) = rwdTimeMatxBin(find(rwdTimeMatxBin < 0)) * -1;

%logistic regression models
glm_rwdTime = fitglm([rwdTimeMatx]', allChoice_R,'distribution','binomial','link','logit'); 
glm_noRwdTime = fitglm([noRwdTimeMatx]', allChoice_R,'distribution','binomial','link','logit');
glm_allTime = fitglm([rwdTimeMatx' noRwdTimeMatx'], allChoice_R,'distribution','binomial','link','logit');


%% create rwds array in time

%find all rwd times in terms of seconds
choiceTimes = ceil(([behSessionData(responseInds).rewardTime] - behSessionData(1).CSon)/1000); %baseline to start time and convert to s from ms
if choiceTimes(1) == 0
    choiceTimes(1) = 1;
end
rwdTimes = choiceTimes(logical(allRewardsBinary));

sessionTime = ceil((behSessionData(end).CSon + 3000 - behSessionData(1).CSon)/1000);     % find total session time and pad time for reward on last trial
sessionRwds = zeros(1,sessionTime);
sessionRwds(rwdTimes) = 1;

%% use exponential decay fit to beta values from time-wise logistic regression to smooth reward hist

xVals = [1:timeBinSize/1000:p.Results.timeMax/1000];
xVals = xVals(1:end-1) + diff(xVals)/2;
expFitTime = singleExpFit(glm_rwdTime.Coefficients.Estimate(2:end), xVals');
expConvTime = expFitTime.a*exp(-(1/expFitTime.b)*(1:p.Results.timeMax/1000));

rwdHxTime = conv(sessionRwds,expConvTime);              %convolve with exponential decay to give weighted moving average
rwdHxTime = rwdHxTime(1:end-(length(expConvTime)-1)); 

choiceTimez = choiceTimes - 1;
if choiceTimez(1) == 0
    rwdHxTimeChoice = [0 rwdHxTime(choiceTimez(2:end))];
else
    rwdHxTimeChoice = rwdHxTime(:,choiceTimez);
end

%%
cumsum_allRchoice = cumsum(allChoices == 1);
cumsum_allLchoice = cumsum(allChoices == -1);
cumsum_allRreward = cumsum(allRewards == 1);
cumsum_allLreward = cumsum(allRewards == -1);
cumsum_blockSwitch = cumsum_allLchoice(blockSwitch);
if cumsum_blockSwitch(1) == 0
    cumsum_blockSwitch(1) = 1;
end

normKern = normpdf(-15:15,0,4);
normKern = normKern / sum(normKern);
halfKern = normKern(round(length(normKern)/2):end);
choiceSlope = atand(diff(conv(cumsum_allRchoice,halfKern))./diff(conv(cumsum_allLchoice,halfKern)));
rwdSlope = atand(diff(conv(cumsum_allRreward,halfKern))./diff(conv(cumsum_allLreward,halfKern)));

avgRwdSlope = [];
tempMax = 0;
for i = 1:length(cumsum_blockSwitch)
    if i ~= length(cumsum_blockSwitch)
        avgRwdSlope(i) = tand(nanmean(rwdSlope(blockSwitch(i):blockSwitch(i+1))));
    else
        avgRwdSlope(i) = tand(mean(rwdSlope(blockSwitch(i):end)));
    end
end

%% make output struct 
s = struct;

s.allChoices = allChoices;
s.allRewards = allRewards;
s.allRewardsBinary = allRewardsBinary;
s.behSessionData = behSessionData;
s.blockSwitch = blockSwitch;
s.choiceHx = choiceHx;
s.choiceTimes = choiceTimes;
s.choiceSlope = choiceSlope;
s.corrLickRate = corrLickRate;
s.lickLat = lickLat;
s.lickLatInds = lickLatInds;
s.lickSide = lickSide;
s.responseInds = responseInds;
s.responseLat = responseLat;
s.responseRateInds = responseRateInds; 
s.rwdHx = rwdHx;
s.rwdHxTimeChoice = rwdHxTimeChoice;
s.rwdSlope = rwdSlope;
s.rwdMatx = rwdMatx;
s.rwdTimeMatx = rwdTimeMatx;
s.rwdTimeMatxBin = rwdTimeMatxBin;
s.sessionRwds = sessionRwds;
s.timeBinSize = timeBinSize;

if p.Results.revForFlag
    s.blockProbs = blockProbs;
end
if ~p.Results.revForFlag
    s.blockSwitchL = blockSwitchL;
    s.blockSwitchR = blockSwitchR;
end


