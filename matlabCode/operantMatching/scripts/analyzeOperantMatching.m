% function analyzeOperantMatching

sessionToImport = '/Volumes/bbari1/BB038/mBB038d20160816/behavior/mBB038d20160816.asc';

sessionText = importMatchingData(sessionToImport);
clear sessionData
sessionData.trialType = [];
sessionData.trialEnd = [];
sessionData.CSon = [];
sessionData.licksL = [];
sessionData.licksR = [];
sessionData.rewardL = [];
sessionData.rewardR = [];
sessionData.rewardTime = [];
sessionData.autoWaterL = [];
sessionData.autoWaterR = [];

for i = 1:length(sessionText)
    % determine beginning and end of trial
    if strfind(sessionText{i},'Trial ') & (regexp(sessionText{i},'Trial ') == 1)% trial begin
        temp1 = regexp(sessionText{i},'('); temp2 = regexp(sessionText{i},')');
        currTrial = str2double(sessionText{i}(temp1(1)+1:temp2(1)-1));
        
        tBegin = i;
        tEndFlag = false;
        j = i + 1;
        while (~tEndFlag)
            if strfind(sessionText{j},'Trial ') & (regexp(sessionText{j},'Trial ') == 1);
                tEnd = j;
                tEndFlag = true;
            else
                j = j + 1;
                if j == length(sessionText)
                    tEnd = length(sessionText);
                    tEndFlag = true;
                end
            end
        end
        
        waterDeliverFlag = false;
        allL_licks = [];
        allR_licks = [];
        for currTrialInd = tBegin+1:tEnd-1
            if strfind(sessionText{currTrialInd},'CS PLUS')
                sessionData(currTrial).trialType = 'CSplus';
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                sessionData(currTrial).CSon = str2double(temp{1}{2});
            elseif strfind(sessionText{currTrialInd},'CS MINUS')
                sessionData(currTrial).trialType = 'CSminus';
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                sessionData(currTrial).CSon = str2double(temp{1}{2});
            end
            if strfind(sessionText{currTrialInd},'L: ') & (regexp(sessionText{currTrialInd},'L: ') == 1);
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                allL_licks = [allL_licks str2double(temp{1}{2})];
            elseif strfind(sessionText{currTrialInd},'R: ') & (regexp(sessionText{currTrialInd},'R: ') == 1);
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                allR_licks = [allR_licks str2double(temp{1}{2})];
            end
            if (~waterDeliverFlag)
                if strfind(sessionText{currTrialInd},'WATER L DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    sessionData(currTrial).rewardL = 1;
                    sessionData(currTrial).rewardR = NaN;
                    sessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                elseif strfind(sessionText{currTrialInd},'WATER L NOT DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    sessionData(currTrial).rewardL = 0;
                    sessionData(currTrial).rewardR = NaN;
                    sessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                elseif strfind(sessionText{currTrialInd},'WATER R DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    sessionData(currTrial).rewardR = 1;
                    sessionData(currTrial).rewardL = NaN;
                    sessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                elseif strfind(sessionText{currTrialInd},'WATER R NOT DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    sessionData(currTrial).rewardR = 0;
                    sessionData(currTrial).rewardL = NaN;
                    sessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                end
            end
            if strfind(sessionText{currTrialInd},'AUTOMATIC WATER L DELIVERED')
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                sessionData(currTrial).autoWaterL = str2double(temp{1}{2});
            elseif strfind(sessionText{currTrialInd},'AUTOMATIC WATER R DELIVERED')
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                sessionData(currTrial).autoWaterR = str2double(temp{1}{2});
            end
            
            if currTrialInd == tEnd-1 || currTrialInd == length(sessionText)-1
                sessionData(currTrial).licksL = allL_licks;
                sessionData(currTrial).licksR = allR_licks;
                if ~waterDeliverFlag
                    sessionData(currTrial).rewardL = NaN;
                    sessionData(currTrial).rewardR = NaN;
                    sessionData(currTrial).rewardTime = NaN;
                end
                if tEnd ~= length(sessionText)
                    temp = regexp(sessionText(tEnd), ': ', 'split');
                    sessionData(currTrial).trialEnd = str2double(temp{1}{2});
                else
                    sessionData(currTrial).trialEnd = NaN;
                end
            end
        end
    end
end

% stringList = {'Trial ', ...
%               'No lick period', ...
%               'CS PLUS', ...
%               'CS MINUS', ...
%               'Odor on', ...
%               'Odor off', ...
%               'ITI - ', ...
%               'WATER L DELIVERED', ...
%               'WATER L NOT DELIVERED', ...
%               'WATER R DELIVERED', ...
%               'WATER R NOT DELIVERED', ...
%               'L:', ...
%               'R:'};
% 
%  outputs = [];
% for currStr = 3:length(stringList)
%     tempInd = strfind(sessionText, stringList{currStr});
%     strInds = find(not(cellfun('isempty', tempInd)));
%     tempExtract = regexp(sessionText(strInds), ': ', 'split');
%     allStr = [tempExtract{:}];
%     outputs{currStr} = cellfun(@str2num,allStr(2:2:end));
% end
% 
% % R = 1, L = -1
% allLicks = sort([outputs{8} outputs{9} outputs{10} outputs{11}]);
% allLickDecisions = allLicks;
% allLickDecisions(ismember(allLickDecisions, outputs{10})) = 1;  % all R licks = 1
% allLickDecisions(ismember(allLickDecisions, outputs{11})) = 1;
% allLickDecisions(ismember(allLickDecisions, outputs{8})) = -1;  % all L licks = -1
% allLickDecisions(ismember(allLickDecisions, outputs{9})) = -1;
% 
% allRewards = allLicks;
% allRewards(ismember(allRewards, outputs{10})) = 1; % all R rewards = 1
% allRewards(ismember(allRewards, outputs{8})) = -1; % all L rewards = 0
% allRewards(allRewards ~=1 & allRewards ~=-1) = 0;
% 
% figure; subplot(3,1,1); hold on; 
% plot(smooth(allLickDecisions(tStart:end),10),'k','linewidth',2)
% plot(smooth(allRewards(tStart:end),10),'b','linewidth',2)
% ylim([-1 1])
% 
% rwdMatx = [];
% for i = 1:10
%     rwdMatx(i,:) = [NaN(1,i) allRewards(tStart:end-i)];
% end
% 
% lickMatx = [];
% for i = 1:10
%     lickMatx(i,:) = [NaN(1,i) allLickDecisions(tStart:end-i)];
% end
% lm = fitlm([rwdMatx]', allLickDecisions(tStart:end),'linear');
% 
% subplot(3,1,2); plot(lm.Coefficients.Estimate(2:end))
% ylim([-0.1 0.8])
% 
% choices = [];
% rewards = [];
% chunkSize = 30;
% for choiceChunk = 1:chunkSize:length(allLickDecisions)-chunkSize
%     choices(((chunkSize-1)+choiceChunk)/chunkSize) = sum(allLickDecisions(choiceChunk:choiceChunk+chunkSize-1) == 1)/chunkSize;
%     rewards(((chunkSize-1)+choiceChunk)/chunkSize) = sum(allRewards(choiceChunk:choiceChunk+chunkSize-1) == 1)/ ...
%         sum(allRewards(choiceChunk:choiceChunk+chunkSize-1) == 1 | allRewards(choiceChunk:choiceChunk+chunkSize-1) == -1);
% end
% 
% subplot(3,1,3); plot(rewards, choices, 'o'); xlabel('Rewards (R > L)'); ylabel('Choices (R > L)');
% xlim([0 1]); ylim([0 1])