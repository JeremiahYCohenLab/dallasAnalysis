function [fisherTbl] = compareChoiceByITI_opMD(xlFile, animals, pre, post, ITIthresh, revForFlag)

if nargin < 5
 ITIthresh = 10;
end

if nargin < 6
 revForFlag = 0;
end

l = ceil(length(animals))/2;
figure;
suptitle(sprintf('win-stay lose-shift by ITI length. thresh = %is', ITIthresh))
for i = 1:length(animals)
    [behTbl_pre, noRwdSwitchTbl_pre, rwdStayTbl_pre] = combineChoiceByITI_opMD(xlFile, animals{i}, pre, ITIthresh, revForFlag);
    [behTbl_post, noRwdSwitchTbl_post, rwdStayTbl_post] = combineChoiceByITI_opMD(xlFile, animals{i}, post, ITIthresh, revForFlag);

    probSwitchPre_s(i) = behTbl_pre.Prob_Switch(1);
    probSwitchPost_s(i) = behTbl_post.Prob_Switch(1);
    probStayPre_s(i) = behTbl_pre.Prob_Stay(1);
    probStayPost_s(i) = behTbl_post.Prob_Stay(1);
    probSwitchPre_l(i) = behTbl_pre.Prob_Switch(2);
    probSwitchPost_l(i) = behTbl_post.Prob_Switch(2);
    probStayPre_l(i) = behTbl_pre.Prob_Stay(2);
    probStayPost_l(i) = behTbl_post.Prob_Stay(2);
    
    noRwdSwitchShort_compare = [noRwdSwitchTbl_pre(1,:); noRwdSwitchTbl_post(1,:)];
    [~, noRwdSwitchShort_p(i), ~] = fishertest(noRwdSwitchShort_compare);
    noRwdSwitchLong_compare = [noRwdSwitchTbl_pre(2,:); noRwdSwitchTbl_post(2,:)]; 
    [~, noRwdSwitchLong_p(i), ~] = fishertest(noRwdSwitchLong_compare);
    rwdStayShort_compare = [rwdStayTbl_pre(1,:); rwdStayTbl_post(1,:)]; 
    [~, rwdStayShort_p(i), ~] = fishertest(rwdStayShort_compare);
    rwdStayLong_compare = [rwdStayTbl_pre(2,:); rwdStayTbl_post(2,:)]; 
    [~, rwdStayLong_p(i), ~] = fishertest(rwdStayLong_compare);
    
%     subplot(2,l,i); hold on;
%     plot([1 2], [behTbl_pre.Prob_Switch(1) behTbl_post.Prob_Switch(1)], '-', 'Color', [1 0.7 1], 'linewidth', 2);
%     plot([1 2], [behTbl_pre.Prob_Switch(2) behTbl_post.Prob_Switch(2)], '-m', 'linewidth', 2);
%     plot([1 2], [behTbl_pre.Prob_Stay(1) behTbl_post.Prob_Stay(1)], '-', 'Color', [0.7 1 1], 'linewidth', 2);
%     plot([1 2], [behTbl_pre.Prob_Stay(2) behTbl_post.Prob_Stay(2)], '-c', 'linewidth', 2);
%     xticks([1 2]); xticklabels({[pre], [post]})
%     xlim([0.75 2.25])
%     if i == 1
%         legend('P(switch | no rwd) short', 'P(switch | no rwd) long', 'P(stay | rwd) short', 'P(stay | rwd) long')
%     end
%     title(animals{i})
    
end

fisherTbl = table(noRwdSwitchShort_p', noRwdSwitchLong_p', rwdStayShort_p', rwdStayLong_p', 'VariableNames',...
    {'noRwdSwitchShort_p', 'noRwdSwitchLong_p', 'rwdStayShort_p', 'rwdStayLong_p'}, 'RowNames', [animals]);

subplot(2,2,1); hold on;
for i = 1:length(animals)
    plot([1 2], [probSwitchPre_s(i) probSwitchPost_s(i)], '-', 'linewidth', 2)
end
xticks([1 2]); xticklabels({[pre], [post]})
xlim([0.75 2.25])
title('P(switch | no rwd) short ITI')
legend([animals])


subplot(2,2,2); hold on;
for i = 1:length(animals)
    plot([1 2], [probSwitchPre_l(i) probSwitchPost_l(i)], '-', 'linewidth', 2)
end
xticks([1 2]); xticklabels({[pre], [post]})
xlim([0.75 2.25])
title('P(switch | no rwd) long ITI')

subplot(2,2,3); hold on;
for i = 1:length(animals)
    plot([1 2], [probStayPre_s(i) probStayPost_s(i)], '-', 'linewidth', 2)
end
xticks([1 2]); xticklabels({[pre], [post]})
xlim([0.75 2.25])
title('P(stay | rwd) short ITI')


subplot(2,2,4); hold on;
for i = 1:length(animals)
    plot([1 2], [probStayPre_l(i) probStayPost_l(i)], '-', 'linewidth', 2)
end
xticks([1 2]); xticklabels({[pre], [post]})
xlim([0.75 2.25])
title('P(stay | rwd) long ITI')
    

 
