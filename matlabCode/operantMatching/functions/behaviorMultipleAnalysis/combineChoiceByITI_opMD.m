 function [behTbl, noRwdSwitchTbl, rwdStayTbl] = combineChoiceByITI_opMD(xlFile, animal, category, ITIthresh, revForFlag)

if nargin < 4
 ITIthresh = 10;
end

if nargin < 5
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


ITIcomb = [];
changeChoiceComb = [];
allRewardsComb = [];

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

    if revForFlag
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
    
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;

    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;
    allRewards = allRewards(1:end-1);
    
    ITI = nan(1,length(responseInds) - 1); 
    for j = 2:length(responseInds)
        ITI(j-1) = behSessionData(responseInds(j)).CSon - behSessionData(responseInds(j-1)).rewardTime;
    end
    ITIcomb = [ITIcomb ITI/1000];
        
    %combine arrays
    changeChoice = [abs(diff(allChoices)) > 0];
    changeChoiceComb = [changeChoiceComb changeChoice];
    allRewardsComb = [allRewardsComb allRewards];  
    
end

shortInds = find(ITIcomb < ITIthresh);
longInds = find(ITIcomb > ITIthresh);
changeChoiceComb_short = changeChoiceComb(shortInds);
stayChoiceComb_short = logical(changeChoiceComb_short == 0);
allRewardsComb_short = allRewardsComb(shortInds);
changeChoiceComb_long = changeChoiceComb(longInds);
stayChoiceComb_long = logical(changeChoiceComb_long == 0);
allRewardsComb_long = allRewardsComb(longInds);

probSwitchNoRwdComb_short = sum(changeChoiceComb_short(allRewardsComb_short==0))/sum(allRewardsComb_short==0);
probSwitchNoRwdComb_long = sum(changeChoiceComb_long(allRewardsComb_long==0))/sum(allRewardsComb_long==0);
probStayRwdComb_short = 1 - (sum(changeChoiceComb_short(allRewardsComb_short==1))/sum(allRewardsComb_short==1));
probStayRwdComb_long = 1 - (sum(changeChoiceComb_long(allRewardsComb_long==1))/sum(allRewardsComb_long==1));

noRwdSwitchTbl = [sum(changeChoiceComb_short(allRewardsComb_short==0)) sum(stayChoiceComb_short(allRewardsComb_short==0));...
    sum(changeChoiceComb_long(allRewardsComb_long==0)) sum(stayChoiceComb_long(allRewardsComb_long==0))];
rwdStayTbl = [sum(stayChoiceComb_short(allRewardsComb_short==1)) sum(changeChoiceComb_short(allRewardsComb_short==1));...
    sum(stayChoiceComb_long(allRewardsComb_long==1)) sum(changeChoiceComb_long(allRewardsComb_long==1))];

behTbl = table([{'short'}; {'long'}], [probSwitchNoRwdComb_short; probSwitchNoRwdComb_long], [probStayRwdComb_short; probStayRwdComb_long], ...
    'VariableNames', {'ITI_duration' 'Prob_Switch' 'Prob_Stay'});


