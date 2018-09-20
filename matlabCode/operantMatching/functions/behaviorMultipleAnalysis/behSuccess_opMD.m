 function [rewardRate, correctRate, noRwdAvg, behTbl] = behSuccess_opMD(xlFile, animal, category, revForFlag)

if nargin < 4
    revForFlag = 0;
end
 
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
normKern = normKern / sum(normKern);

choiceRange = 15;
noRwdCombined = [];
changeHistogram = [];
rwdHistChange = []; 
changeChoiceComb = [];
allRewardsComb = [];
allRewardsBinComb = [];

for i = 1: length(dayList)
    sessionName = dayList{i};
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
        if revForFlag
            behSessionData = sessionData;
        end
    elseif revForFlag                                    %otherwise generate the struct
        [behSessionData, ~] = generateSessionData_behav_operantMatching(sessionName);
    else
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
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
    x = [ones(1, tmpCorrect) zeros(1, (tmpResp - tmpCorrect))];
    y = [ones(1, tmpResp/2) zeros(1, tmpResp/2)];
    randChance(i) = ranksum(x,y);
    
    %find fraction of trials that are rewarded and fraction that are on the higher spout
    allRewardsBin = allRewards;
    allRewardsBin(allRewards == -1) = 1;
    rewardFraction(i) = sum(allRewardsBin)/length(allRewardsBin);
    correctFraction(i) = correctFraction(i)/length(allRewards);
    allRewardsBinComb = [allRewardsBinComb allRewardsBin];
    
    %find distance between smoothed rewards and choices   
    smoothChoices = conv(allChoices,normKern, 'same')/max(conv(allChoices,normKern,'same'));
    smoothRewards = conv(allRewards,normKern,'same')/max(conv(allRewards,normKern,'same'));
    avgDist(i) = mean(abs(smoothChoices - smoothRewards));
    semDist(i) = std(abs(smoothChoices - smoothRewards)) / sqrt(length(allRewards));
    

    %find choice behavior prior to a zero crossing of smoothed choices
%     LtoR = [];
%     RtoL = [];
%     for j = 1:length(smoothChoices)-1
%         if smoothChoices(j) < 0 & smoothChoices(j+1) > 0
%             LtoR = [LtoR j];
%         elseif smoothChoices(j) > 0 & smoothChoices(j+1) < 0
%             RtoL = [RtoL j];
%         end
%     end
%     for j = 1:length(LtoR)
%         if LtoR(j) > choiceRange
%             noRwdL(j) = length(find(allRewards(LtoR(j)-choiceRange:LtoR(j)) == 0 & allChoices(LtoR(j)-choiceRange:LtoR(j)) == -1));
%         else
%             noRwdL(j) = NaN;
%         end
%     end
%     for j = 1:length(RtoL)
%         if RtoL(j) > choiceRange
%             noRwdR(j) = length(find(allRewards(RtoL(j)-choiceRange:RtoL(j)) == 0 & allChoices(RtoL(j)-choiceRange:RtoL(j)) == 1));
%         else
%             noRwdR(j) = NaN;
%         end
%     end
%     avgNoRwd(i) = nanmean([noRwdR noRwdL]);
%     noRwdCombined = [noRwdCombined noRwdR noRwdL];
        
    %find probability of switch/stay given no rwd/rwd
    changeChoice = [abs(diff(allChoices)) > 0];
    allRewardz = allRewardsBin(1:end-1);
    probSwitchNoRwd(i) = sum(changeChoice(allRewardz==0))/sum(allRewardz==0);
    probStayRwd(i) = 1 - (sum(changeChoice(allRewardz==1))/sum(allRewardz==1));
    normSwitches(i) = sum(changeChoice)/length(allChoices);
    
    changeChoiceComb = [changeChoiceComb changeChoice];
    allRewardsComb = [allRewardsComb allRewardz];
 
   %% 
%     tMax = 10;
%     rwdMatx = [];
%     for j = 1:tMax
%         rwdMatx(j,:) = [NaN(1,j) allRewards(1:end-j)];
%     end
% 
%     glm_rwd = fitglm([rwdMatx]', allChoice_R,'distribution','binomial','link','logit');
%     
%     expFit = singleExpFit(glm_rwd.Coefficients.Estimate(2:end));
%     expConv = expFit.a*exp(-(1/expFit.b)*(1:10));
%     expConv = expConv./sum(expConv);
% 
%     rwdsTemp = find(allRewards == -1);                          %make all rewards have the same value
%     allRewardsBinary = allRewards;
%     allRewardsBinary(rwdsTemp) = 1;
%     rwdHx = conv(allRewardsBinary,expConv);              %convolve with exponential decay to give weighted moving average
%     rwdHx = rwdHx(1:end-(length(expConv)-1));                   %to account for convolution padding
%     rwdHx_L = conv(allReward_L,expConv);                 %same convolution but only with L rewards over all trials
%     rwdHx_L = rwdHx_L(1:end-(length(expConv)-1));  
%     rwdHx_R = conv(allReward_R,expConv);
%     rwdHx_R = rwdHx_R(1:end-(length(expConv)-1));
%     
%     for j = find(changeChoice == 1)
%         if j > 1
%             if allChoices(j) == 1 
%                 temp = 0;
%                 goBack = 1;
%                 while (j - goBack > 0) && allChoices(j-goBack) == -1 && allRewards(j-goBack) == 0
%                     temp = temp + 1;
%                     goBack = goBack + 1;
%                 end
%                 changeHistogram = [changeHistogram temp];
%                 rwdHistChange = [rwdHistChange rwdHx(j-1)];
%             elseif allChoices(j) == -1
%                 temp = 0;
%                 goBack = 1;
%                 while (j - goBack > 0) && allChoices(j-goBack) == 1 && allRewards(j-goBack) == 0 % 
%                     temp = temp + 1;
%                     goBack = goBack + 1;
%                 end
%                 changeHistogram = [changeHistogram temp];
%                 rwdHistChange = [rwdHistChange rwdHx(j-1)];
%             end
%         end
%     end
    
    
%%
    ITIlicks(i) = 0;
    for j = 1:length(behSessionData)
        ITIlicks(i) =  ITIlicks(i) + sum([behSessionData(j).licksR > (behSessionData(j).CSon + 2500) behSessionData(j).licksL > (behSessionData(j).CSon + 2500)]);
    end
    ITIlicks(i) = ITIlicks(i)/length(behSessionData);
    
%%
    %find model parameters
%     x = qLearning_fit2LR([dayList{i} '.asc']);
%     alphaNPE(i) = x.threeParams_twoLearnRates.bestParams(1);
%     alphaPPE(i) = x.threeParams_twoLearnRates.bestParams(2);
%     beta(i) = x.threeParams_twoLearnRates.bestParams(3);
    
end

probSwitchNoRwdComb = sum(changeChoiceComb(allRewardsComb==0))/sum(allRewardsComb==0);
probStayRwdComb = 1 - (sum(changeChoiceComb(allRewardsComb==1))/sum(allRewardsComb==1));

correctRate = correct/(correct + incorrect);
rewardRate = sum(allRewardsBinComb)/length(allRewardsBinComb);
noRwdAvg = nanmean(noRwdCombined);

% behTbl = table(dayList, correctFraction', rewardFraction', avgDist', semDist', probSwitchNoRwd', probStayRwd', randChance', normSwitches', avgNoRwd', ITIlicks',...
%     'VariableNames', {'Session' 'Fraction_Correct' 'Fraction_Rewarded' 'Distance_Avg' 'Distance_SEM' 'Prob_Switch' 'Prob_Stay' 'pVal_Chance'... 
%     'Norm_Switches' 'No_Reward' 'ITI_Licks'});

behTbl = table(dayList, correctFraction', rewardFraction', avgDist', semDist', probSwitchNoRwd', probStayRwd', randChance', normSwitches', ITIlicks',...
    'VariableNames', {'Session' 'Fraction_Correct' 'Fraction_Rewarded' 'Distance_Avg' 'Distance_SEM' 'Prob_Switch' 'Prob_Stay' 'pVal_Chance'... 
    'Norm_Switches' 'ITI_Licks'});

