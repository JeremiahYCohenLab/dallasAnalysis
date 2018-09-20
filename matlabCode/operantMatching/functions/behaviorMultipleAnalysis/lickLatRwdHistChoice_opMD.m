 function [behTbl, mdl] = lickLatRwdHistChoice_opMD(xlFile, animal, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag',0)
p.addParameter('tMax', 20)
p.addParameter('timeMax', 121000)
p.addParameter('timeBins', 12)
p.parse(varargin{:}); 

timeBinSize = (p.Results.timeMax - 1000)/p.Results.timeBins;
timeBinEdges = [1000:timeBinSize:p.Results.timeMax];
 
[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
correct = 0;
incorrect = 0;
correctFraction = zeros(1, length(dayList));
rewardFraction = zeros(1, length(dayList));

%guassian kernel for smoothing of raw rewards and choices
normKern = normpdf(-5:5,0,4);
normKern = normKern / sum(normKern);                                %make sure this produces a causal filter (i think need to split this in half)


for i = 1: length(dayList)
    fprintf('Session number: %d of %d \n', i, length(dayList));
    sessionName = dayList{i};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName(end) sep sessionName '_sessionData.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sessionName '_sessionData.mat'];
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
    
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    allChoice_R = double(allChoices == 1);
    allChoice_L = double(allChoices == -1);
    
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;

    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = -1;
    for j = 1:length(responseInds)
        if (allChoices(j) == 1 & (behSessionData(responseInds(j)).rewardProbR >= behSessionData(responseInds(j)).rewardProbL) | ...
                allChoices(j) == -1 & (behSessionData(responseInds(j)).rewardProbL >= behSessionData(responseInds(j)).rewardProbR))
            correct = correct + 1;
            correctFraction(i) = correctFraction(i) + 1;
        else
            incorrect = incorrect + 1;
        end    
    end
    
    %see if behavior differs from chance
    if mod(length(responseInds),2) == 1
        tmpCorrect = correctFraction(i) - 1;
        tmpResp = length(responseInds) - 1;
    else
        tmpCorrect = correctFraction(i);
        tmpResp = length(responseInds);
    end
    mdl = [ones(1, tmpCorrect) zeros(1, (tmpResp - tmpCorrect))];
    y = [ones(1, tmpResp/2) zeros(1, tmpResp/2)];
    randChance(i) = ranksum(mdl,y);
    
    %find fraction of trials that are rewarded and fraction that are on the higher spout
    allRewardsBin = allRewards;
    allRewardsBin(allRewards == -1) = 1;
    rewardFraction(i) = sum(allRewardsBin)/length(allRewardsBin);
    correctFraction(i) = correctFraction(i)/length(allRewards);
    
    %find distance between smoothed rewards and choices   
    smoothChoices = conv(allChoices,normKern, 'same')/max(conv(allChoices,normKern,'same'));
    smoothRewards = conv(allRewards,normKern,'same')/max(conv(allRewards,normKern,'same'));
    avgDist(i) = mean(abs(smoothChoices - smoothRewards));
    semDist(i) = std(abs(smoothChoices - smoothRewards)) / sqrt(length(allRewards));
    

    %find choice behavior prior to a zero crossing of smoothed choices
    LtoR = [];
    RtoL = [];
    for j = 1:length(smoothChoices)-1
        if smoothChoices(j) < 0 & smoothChoices(j+1) > 0
            LtoR = [LtoR j];
        elseif smoothChoices(j) > 0 & smoothChoices(j+1) < 0
            RtoL = [RtoL j];
        end
    end
    for j = 1:length(LtoR)
        if LtoR(j) > 10
            noRwdL(j) = length(find(allRewards(LtoR(j)-10:LtoR(j)) == 0 & allChoices(LtoR(j)-10:LtoR(j)) == -1));
        else
            noRwdL(j) = NaN;
        end
    end
    for j = 1:length(RtoL)
        if RtoL(j) > 10
            noRwdR(j) = length(find(allRewards(RtoL(j)-10:RtoL(j)) == 0 & allChoices(RtoL(j)-10:RtoL(j)) == 1));
        else
            noRwdR(j) = NaN;
        end
    end
    avgNoRwd(i) = nanmean([noRwdR noRwdL]);
        
    %find probability of switch/stay given no rwd/rwd
    changeChoice = [abs(diff(allChoices)) > 0];
    allRewards(allRewards==-1) = 1;
    allRewards = allRewards(1:end-1);
    probSwitchNoRwd(i) = sum(changeChoice(allRewards==0))/sum(allRewards==0);
    probStayRwd(i) = 1 - (sum(changeChoice(allRewards==1))/sum(allRewards==1));
    normSwitches(i) = sum(changeChoice)/length(allChoices);
    
    %find model parameters
    mdl = qLearning_fit2LR([dayList{i} '.asc']);
    alphaNPE(i) = mdl.fiveParams_twoLearnRates_alphaForget_bias.bestParams(1);
    alphaPPE(i) = mdl.fiveParams_twoLearnRates_alphaForget_bias.bestParams(2);
    alphaForget(i) = mdl.fiveParams_twoLearnRates_alphaForget_bias.bestParams(3);
    beta(i) = mdl.fiveParams_twoLearnRates_alphaForget_bias.bestParams(4);
    bias(i) = mdl.fiveParams_twoLearnRates_alphaForget_bias.bestParams(5);
    
    
    %% compare rewards in time to lick latencies
    
    rwdMatx = zeros(p.Results.tMax, length(responseInds)); 
    for j = 2:length(responseInds)          
        k = 1;
        %find time between "current" choice and previous rewards, up to p.Results.timeMax in the past 
        timeTmpL = []; timeTmpR = [];
        while j-k > 0 & behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime < p.Results.timeMax
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
            for k = 1:p.Results.tMax
                if ~isempty(binnedRwds == k)
                    rwdMatx(k,j) = -1*sum(binnedRwds == k);
                else
                    rwdMatx(k,j) = NaN;
                end
            end
        end
        if ~isempty(timeTmpR)
            binnedRwds = discretize(timeTmpR,timeBinEdges);
            for k = 1:p.Results.tMax
                if ~isempty(binnedRwds == k) & isnan(rwdMatx(k,j))
                    rwdMatx(k,j) = sum(binnedRwds == k);
                elseif ~isempty(binnedRwds == k) & ~isnan(rwdMatx(k,j))
                    rwdMatx(k,j) = rwdMatx(k,j) + sum(binnedRwds == k);
                else
                    rwdMatx(k,j) = NaN;
                end
            end
        end
    end
    
    %fill in NaNs at beginning of session
    j = 2;
    while behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(1)).rewardTime < p.Results.timeMax
        tmpDiff = behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(1)).rewardTime;
        binnedDiff = discretize(tmpDiff, timeBinEdges);
        rwdMatx(binnedDiff:p.Results.tMax,j) = NaN;
        j = j+1;
    end
    rwdMatx(:,1) = NaN;
    rwdBinMatx = rwdMatx; rwdBinMatx(find(rwdBinMatx == -1)) = 1;
    
    lickLat = [behSessionData(responseInds).rewardTime] - [behSessionData(responseInds).CSon];
    indsR = find(allChoices == 1);
    indsL = find(allChoices == -1);
    lickLat_R = zscore(lickLat(indsR));
    lickLat_L = zscore(lickLat(indsL));
    lickLat = NaN(1, length(allChoices));
    lickLat(indsR) = lickLat_R;
    lickLat(indsL) = lickLat_L;
    
    lickLat = [NaN lickLat(2:end)];

    %linear regression model to find time constant of reward effect
    lm_rwdLick = fitlm([rwdBinMatx]', lickLat);
    expFit = singleExpFit(lm_rwdLick.Coefficients.Estimate(2:end));
    tau(i) = 1/expFit.b;
    
    %linear regression to find relationship between choice-relevant reward hist and lick latency
    glm_rwdChoice = fitglm([rwdMatx]', allChoice_R, 'distribution','binomial','link','logit');
    expFit = singleExpFit(glm_rwdChoice.Coefficients.Estimate(2:end));
    expConv = expFit.a*exp(-(1/expFit.b)*(1:10));
    expConv = expConv./sum(expConv);
    rwdHx = conv(allRewardsBin,expConv);              %convolve with exponential decay to give weighted moving average
    rwdHx = rwdHx(1:end-(length(expConv)-1));            %to account for convolution padding
    rwdHx = [0 rwdHx(1:end-1)];                          %to make rwdHx for each trial prior to its outcome
    
    lm_rwdHxLick = fitlm(rwdHx', lickLat);
    rSqr(i) = lm_rwdHxLick.Rsquared.Adjusted;
    

end

correctRate = correct/(correct + incorrect);

behTbl = table(dayList, correctFraction', rewardFraction', avgDist', semDist', probSwitchNoRwd', probStayRwd', randChance',...
    normSwitches', avgNoRwd', tau', rSqr', alphaNPE', alphaPPE', alphaForget', beta', bias',...
    'VariableNames', {'Session' 'Fraction_Correct' 'Fraction_Rewarded' 'Distance_Avg' 'Distance_SEM' 'Prob_Switch' 'Prob_Stay'...
    'pVal_Chance' 'Norm_Switches' 'No_Reward' 'Rwd_Lick_Tau', 'Choice_Rwd_Lick_R' 'Alpha_NPE' 'Alpha_PPE' 'Alpha_Forget'...
    'Beta' 'Bias'});

mdl = fitlm(lickTbl_CG14_effect,'Choice_Rwd_Lick_R ~ Fraction_Correct + Fraction_Rewarded + Distance_Avg + Distance_SEM + Prob_Switch + Prob_Stay + pVal_Chance + Norm_Switches + No_Reward + Rwd_Lick_Tau + Alpha_NPE + Alpha_PPE + Alpha_Forget + Beta + Bias')