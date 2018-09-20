function [fourParams, oppoParams, rBar] = QlearningModelComparisonSession_opMD(xlFile, animal, category, revForFlag)

if nargin < 4
    revForFlag = 0;
end

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

rBar = [];
%loop for each session in the list
for i = 1: length(dayList)
    fprintf('File number: %d of %d \n', i, length(dayList));
    sessionName = dayList{i};                       %extract relevant info from session title
    
    Qmdls = qLearning_fitOpponency([sessionName '.asc'], revForFlag);
    
    fourParams(i,:) = [Qmdls.fourParams_twoLearnRates_alphaForget.bestParams Qmdls.fourParams_twoLearnRates_alphaForget.BIC];
    oppoParams(i,:) = [Qmdls.fiveParams_opponency.bestParams Qmdls.fiveParams_opponency.BIC];
    rBar = [rBar Qmdls.fiveParams_opponency.rBar'];
    
end

% figure; suptitle('4 param')
% subplot(1,4,1); hold on; title('alpha NPE')
% errorbar([1 2], [mean(CG33pre_fourParams(:,1)) mean(CG33post_fourParams(:,1))],...
%     [std(CG33pre_fourParams(:,1))/sqrt(length(CG33pre_fourParams(:,1)))  std(CG33post_fourParams(:,1))/sqrt(length(CG33post_fourParams(:,1)))],'-c', 'linewidth', 2)
% errorbar([1 2], [mean(CG34pre_fourParams(:,1)) mean(CG34post_fourParams(:,1))],...
%     [std(CG34pre_fourParams(:,1))/sqrt(length(CG34pre_fourParams(:,1)))  std(CG34post_fourParams(:,1))/sqrt(length(CG34post_fourParams(:,1)))],'-m', 'linewidth', 2)
% xlim([0.75 2.25])
% legend('CG33', 'CG34')
% 
% subplot(1,4,2); hold on; title('alpha PPE')
% errorbar([1 2], [mean(CG33pre_fourParams(:,2)) mean(CG33post_fourParams(:,2))],...
%     [std(CG33pre_fourParams(:,2))/sqrt(length(CG33pre_fourParams(:,2)))  std(CG33post_fourParams(:,2))/sqrt(length(CG33post_fourParams(:,2)))],'-c', 'linewidth', 2)
% errorbar([1 2], [mean(CG34pre_fourParams(:,2)) mean(CG34post_fourParams(:,2))],...
%     [std(CG34pre_fourParams(:,2))/sqrt(length(CG34pre_fourParams(:,2)))  std(CG34post_fourParams(:,2))/sqrt(length(CG34post_fourParams(:,2)))],'-m', 'linewidth', 2)
% xlim([0.75 2.25])
% 
% subplot(1,4,3); hold on; title('alpha forget')
% errorbar([1 2], [mean(CG33pre_fourParams(:,3)) mean(CG33post_fourParams(:,3))],...
%     [std(CG33pre_fourParams(:,3))/sqrt(length(CG33pre_fourParams(:,3)))  std(CG33post_fourParams(:,3))/sqrt(length(CG33post_fourParams(:,3)))],'-c', 'linewidth', 2)
% errorbar([1 2], [mean(CG34pre_fourParams(:,3)) mean(CG34post_fourParams(:,3))],...
%     [std(CG34pre_fourParams(:,3))/sqrt(length(CG34pre_fourParams(:,3)))  std(CG34post_fourParams(:,3))/sqrt(length(CG34post_fourParams(:,3)))],'-m', 'linewidth', 2)
% xlim([0.75 2.25])
% 
% subplot(1,4,4); hold on; title('beta')
% errorbar([1 2], [mean(CG33pre_fourParams(:,4)) mean(CG33post_fourParams(:,4))],...
%     [std(CG33pre_fourParams(:,4))/sqrt(length(CG33pre_fourParams(:,4)))  std(CG33post_fourParams(:,4))/sqrt(length(CG33post_fourParams(:,4)))],'-c', 'linewidth', 2)
% errorbar([1 2], [mean(CG34pre_fourParams(:,4)) mean(CG34post_fourParams(:,4))],...
%     [std(CG34pre_fourParams(:,4))/sqrt(length(CG34pre_fourParams(:,4)))  std(CG34post_fourParams(:,4))/sqrt(length(CG34post_fourParams(:,4)))],'-m', 'linewidth', 2)
% xlim([0.75 2.25])

