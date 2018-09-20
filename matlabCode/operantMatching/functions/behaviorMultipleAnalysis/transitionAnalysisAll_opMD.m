function [transChoiceMatxForty, transChoiceMatxSeventy, mdlFitForty, mdlFitSeventy] = transitionAnalysisAll_opMD(xlFile, sheets)

[root, sep] = currComputer();
transChoiceMatxForty = []; 
transChoiceMatxSeventy = [];
range = 25;

for k = 1:length(sheets)
    [~, dayList, ~] = xlsread(xlFile, sheets{k});
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
        else
            [behSessionData, blockSwitch, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
        end

        responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
        omitInds = isnan([behSessionData.rewardTime]); 

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

        allReward_R(isnan(allReward_R)) = 0;
        allReward_L(isnan(allReward_L)) = 0;
        allChoice_R = double(allChoices == 1);
        allChoice_L = double(allChoices == -1);

        allRewards = zeros(1,length(allChoices));
        allRewards(logical(allReward_R)) = 1;
        allRewards(logical(allReward_L)) = -1;

        rwdProb_R = [behSessionData(responseInds).rewardProbR]; 
        rwdProb_L = [behSessionData(responseInds).rewardProbL]; 

        for j = 2:(length(blockSwitch) - 1)
            tmpInd = blockSwitch(j);
            if rwdProb_R(tmpInd) == 70 & rwdProb_L(tmpInd) == 10 & rwdProb_R(tmpInd+1) == 10 & rwdProb_L(tmpInd+1) == 70
    %            if sum(allChoice_R((tmpInd - 9):tmpInd)) > 7
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transChoiceMatxSeventy = [transChoiceMatxSeventy; allChoices((tmpInd-range+1):(tmpInd+range))];
                    end
    %            end
            elseif rwdProb_R(tmpInd) == 40 & rwdProb_L(tmpInd) == 10 & rwdProb_R(tmpInd+1) == 10 & rwdProb_L(tmpInd+1) == 70
    %            if sum(allChoice_R((tmpInd - 9):tmpInd)) > 7
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transChoiceMatxForty = [transChoiceMatxForty; allChoices((tmpInd-range+1):(tmpInd+range))];
                    end
    %            end
            elseif rwdProb_L(tmpInd) == 70 & rwdProb_R(tmpInd) == 10 & rwdProb_L(tmpInd+1) == 10 & rwdProb_R(tmpInd+1) == 70
    %            if sum(allChoice_L((tmpInd - 9):tmpInd)) > 7
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transChoiceMatxSeventy = [transChoiceMatxSeventy; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                    end
    %            end        
            elseif rwdProb_L(tmpInd) == 40 & rwdProb_R(tmpInd) == 10 & rwdProb_L(tmpInd+1) == 10 & rwdProb_R(tmpInd+1) == 70
    %            if sum(allChoice_L((tmpInd - 9):tmpInd)) > 7
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transChoiceMatxForty = [transChoiceMatxForty; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                    end
    %            end
            end   
        end
    end   
end

fortyAvg = mean(transChoiceMatxForty,1);
seventyAvg = mean(transChoiceMatxSeventy,1);

figure; hold on
x = [-range+1:range];
plot(x,fortyAvg, 'r')
plot(x,seventyAvg, 'b')
title([ 'Choice at block transitions to 10/70'])
ylabel('Choice average')
xlabel('Trials from switch')
legend('40 -> 10', '70 -> 10')
linetype = 'k'
vline(0, linetype)

xx = [1:range+1];
mdlFitForty = singleExpFitInt(xx,fortyAvg(range:end));
mdlFitSeventy = singleExpFitInt(xx,seventyAvg(range:end));
expConvForty = mdlFitForty.a*exp(-(mdlFitForty.b)*(1:range+1))+mdlFitForty.c;
expConvSeventy = mdlFitSeventy.a*exp(-(mdlFitSeventy.b)*(1:range+1))+mdlFitSeventy.c;
plot([0:25], expConvForty, '--r'); plot([0:25], expConvSeventy,'--b')
ylim([-1 1])
    