function wslsRwdHx_opMD(xlFile, animals, categories, revForFlag, plotFlag)

if nargin < 5
    plotFlag = 0;
end
if nargin < 4
    revForFlag = 0;
end

[root, sep] = currComputer();
tMax = 12;

for i = 1:length(animals)
    
    [~, dayList, ~] = xlsread(xlFile, animals{i});
    [~,col] = find(~cellfun(@isempty,strfind(dayList, categories{i})) == 1);
    dayList = dayList(2:end,col);
    endInd = find(cellfun(@isempty,dayList),1);
    if ~isempty(endInd)
        dayList = dayList(1:endInd-1,:);
    end
    
    [glm, tMax] = combineLogReg_opMD(xlFile, animals{i}, categories{i}, revForFlag);
    expFit = singleExpFit(glm.Coefficients.Estimate(2:tMax+1));
    expConv = expFit.a*exp(-(1/expFit.b)*(1:tMax));
    expConv = expConv./sum(expConv);

    combinedRwdHx = [];
    combinedRwds = [];
    combinedChangeChoice = [];
    
    for j = 1: length(dayList)
        sessionName = dayList{j};
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
            [behSessionData, ~] = generateSessionData_operantMatching(sessionName);
        else
            [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
        end

        responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
        allReward_R = [behSessionData(responseInds).rewardR]; 
        allReward_L = [behSessionData(responseInds).rewardL]; 
        allChoices = NaN(1,length(behSessionData(responseInds)));
        allChoices(~isnan(allReward_R)) = 1;
        allChoices(~isnan(allReward_L)) = -1;

        allReward_R(isnan(allReward_R)) = 0;
        allReward_L(isnan(allReward_L)) = 0;
        allChoice_R = double(allChoices == 1);
        allChoice_L = double(allChoices == -1);

        allRewards = zeros(1,length(allChoices));
        allRewards(logical(allReward_R)) = 1;
        allRewards(logical(allReward_L)) = -1;
        allRewardsBin = allRewards;
        allRewardsBin(allRewards == -1) = 1;

        allNoRewards = allChoices;
        allNoRewards(logical(allReward_R)) = 0;
        allNoRewards(logical(allReward_L)) = 0;
    
        rwdHx = conv(allRewardsBin, expConv);
        rwdHx = rwdHx(1:(end-(length(expConv)-1)));
        
        combinedRwdHx = [combinedRwdHx rwdHx(1:end-2)];
        combinedRwds = [combinedRwds allRewardsBin(2:end-1)];
        changeChoice = [abs(diff(allChoices)) > 0];
        combinedChangeChoice = [combinedChangeChoice changeChoice(2:end)];
        
    end
    
    rwdHxInds_low = logical(combinedRwdHx < 1/3);
    rwdHxInds_med = logical(combinedRwdHx > 1/3 & combinedRwdHx < 2/3);
    rwdHxInds_high = logical(combinedRwdHx > 2/3);
    probSwitchNoRwd_low(i) = sum(combinedChangeChoice(combinedRwds==0 & rwdHxInds_low))/sum(combinedRwds==0 & rwdHxInds_low);
    probSwitchNoRwd_med(i) = sum(combinedChangeChoice(combinedRwds==0 & rwdHxInds_med))/sum(combinedRwds==0 & rwdHxInds_med);
    probSwitchNoRwd_high(i) = sum(combinedChangeChoice(combinedRwds==0 & rwdHxInds_high))/sum(combinedRwds==0 & rwdHxInds_high);
    probStayRwd_low(i) = 1 - (sum(combinedChangeChoice(combinedRwds==1 & rwdHxInds_low))/sum(combinedRwds==1 & rwdHxInds_low));
    probStayRwd_med(i) = 1 - (sum(combinedChangeChoice(combinedRwds==1 & rwdHxInds_med))/sum(combinedRwds==1 & rwdHxInds_med));
    probStayRwd_high(i) = 1 - (sum(combinedChangeChoice(combinedRwds==1 & rwdHxInds_high))/sum(combinedRwds==1 & rwdHxInds_high));
    

end

figure;
subplot(1,2,1); hold on;
for i = 1:length(animals)
    plot([1 2 3], [probStayRwd_low(i) probStayRwd_med(i) probStayRwd_high(i)], 'LineWidth', 2)
end
xlim([0.5 3.5])
xticklabels({'low', 'medium', 'high'})
ylabel('probability')
title('win-stay')
legend([animals])
subplot(1,2,2); hold on;
for i = 1:length(animals)
    plot([1 2 3], [probSwitchNoRwd_low(i) probSwitchNoRwd_med(i) probSwitchNoRwd_high(i)], 'LineWidth', 2)
end
xlim([0.5 3.5])
xticklabels({'low', 'medium', 'high'})
ylabel('probability')
title('lose-shift')
end