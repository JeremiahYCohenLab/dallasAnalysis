function [aicTrialAvg, aicTimeAvg, semTrial, semTime] = logRegModelComparisonSession_opMD(xlFile, animal, category)

%determine root for file location
[root, sep] = currComputer();

%import behavior session titles for desired category
[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

tMax = 20;
timeMax = 127000;
binSize = (timeMax - 1000)/tMax;
timeBinEdges = [1000:binSize:timeMax];  %no trials shorter than 1s between outcome and CS on

%loop for each session in the list
for i = 1: length(dayList)              
    sessionName = dayList{i};                       %extract relevant info from session title
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')         %define appropriate data path
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')        %load preprocessed struct if there is one
        load(sessionDataPath)
    else                                    %otherwise generate the struct
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
    
    %create arrays for choices and rewards
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
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
    allRewards(logical(allReward_L)) = 1;
    

    rwdMatx = [];
    for j = 1:tMax
        rwdMatx(j,:) = [NaN(1,j) allRewards(1:end-j)];
    end
    glm_rwd = fitglm([rwdMatx]', allChoice_R,'distribution','binomial','link','logit');
    aicTrial(i) = aicbic(glm_rwd.LogLikelihood, tMax);
    
    
    %logistic regression for rwds in time
    rwdTimeMatx = NaN(tMax, length(responseInds));     %initialize matrices for number of response trials x number of time bins
    noRwdTimeMatx =  NaN(tMax, length(responseInds));
    for j = 2:length(responseInds)          
        k = 1;
        %find time between "current" choice and previous rewards, up to timeMax in the past 
        timeTmpL = []; timeTmpR = []; nTimeTmpL = []; nTimeTmpR = [];
        while j-k > 0 & behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime < timeMax
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
            for k = 1:tMax
                if ~isempty(binnedRwds == k)
                    rwdTimeMatx(k,j) = -1*sum(binnedRwds == k);
                else
                    rwdTimeMatx(k,j) = NaN;
                end
            end
        end
        if ~isempty(timeTmpR)
            binnedRwds = discretize(timeTmpR,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedRwds == k) & isnan(rwdTimeMatx(k,j))
                    rwdTimeMatx(k,j) = sum(binnedRwds == k);
                elseif ~isempty(binnedRwds == k) & ~isnan(rwdTimeMatx(k,j))
                    rwdTimeMatx(k,j) = rwdTimeMatx(k,j) + sum(binnedRwds == k);
                else
                    rwdTimeMatx(k,j) = NaN;
                end
            end
        end
        if ~isempty(nTimeTmpL)
            binnedNoRwds = discretize(nTimeTmpL,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedNoRwds == k)
                    noRwdTimeMatx(k,j) = -1*sum(binnedNoRwds == k);
                else
                    noRwdTimeMatx(k,j) = NaN;
                end
            end
        end
        if ~isempty(nTimeTmpR)
            binnedNoRwds = discretize(nTimeTmpR,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedNoRwds == k) & isnan(noRwdTimeMatx(k,j))
                    noRwdTimeMatx(k,j) = sum(binnedNoRwds == k);
                elseif ~isempty(binnedNoRwds == k) & ~isnan(noRwdTimeMatx(k,j))
                    noRwdTimeMatx(k,j) = noRwdTimeMatx(k,j) + sum(binnedNoRwds == k);
                else
                    noRwdTimeMatx(k,j) = NaN;
                end
            end
        end
    end

    rwdTimeMatx(:,1) = NaN;
    glm_rwdTime = fitglm([rwdTimeMatx]', allChoice_R,'distribution','binomial','link','logit');
%    glm_all = fitglm([rwdTimeMatx' noRwdTimeMatx'], allChoice_R,'distribution','binomial','link','logit');
    aicTime(i) = aicbic(glm_rwdTime.LogLikelihood, tMax);
    
end

aicTrialAvg = mean(aicTrial); 
aicTimeAvg = mean(aicTime); 
semTrial = std(aicTrial)/sqrt(length(aicTrial));
semTime = std(aicTime)/sqrt(length(aicTime));

% figure; hold on;
% errorbar([1 2],[aicTrialAvg aicTimeAvg],[semTrial semTime],[semTrial semTime],'Color', [0.7 0 1],'linewidth',2)
% ylabel('mean AIC score')
% xlim([0.5 2.5])
% xticks([1 2])
% xticklabels({'trial', 'time'})
% title([animal ' ' category])

