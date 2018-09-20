figure;

subplot(1,4,1); hold on;
plot([1 2], [mdlCG05pre.fourParams_twoLearnRates_alphaForget.bestParams(1) mdlCG05post.fourParams_twoLearnRates_alphaForget.bestParams(1)],'r','linewidth',2);
plot([1 2], [mdlCG14preS.fourParams_twoLearnRates_alphaForget.bestParams(1) mdlCG14effect.fourParams_twoLearnRates_alphaForget.bestParams(1)],'b','linewidth',2);
plot([1 2], [mdlCG15preS.fourParams_twoLearnRates_alphaForget.bestParams(1) mdlCG15effect.fourParams_twoLearnRates_alphaForget.bestParams(1)],'Color', [0.7 0 1],'linewidth',2);
xlim([0.75 2.25])
xticks([1 2])
xticklabels({'pre','post'})
title('alpha NPE')

subplot(1,4,2); hold on;
plot([1 2], [mdlCG05pre.fourParams_twoLearnRates_alphaForget.bestParams(2) mdlCG05post.fourParams_twoLearnRates_alphaForget.bestParams(2)],'r','linewidth',2);
plot([1 2], [mdlCG14preS.fourParams_twoLearnRates_alphaForget.bestParams(2) mdlCG14effect.fourParams_twoLearnRates_alphaForget.bestParams(2)],'b','linewidth',2);
plot([1 2], [mdlCG15preS.fourParams_twoLearnRates_alphaForget.bestParams(2) mdlCG15effect.fourParams_twoLearnRates_alphaForget.bestParams(2)],'Color', [0.7 0 1],'linewidth',2);
xlim([0.75 2.25])
xticks([1 2])
xticklabels({'pre','post'})
title('alpha PPE')

subplot(1,4,3); hold on;
plot([1 2], [mdlCG05pre.fourParams_twoLearnRates_alphaForget.bestParams(3) mdlCG05post.fourParams_twoLearnRates_alphaForget.bestParams(3)],'r','linewidth',2);
plot([1 2], [mdlCG14preS.fourParams_twoLearnRates_alphaForget.bestParams(3) mdlCG14effect.fourParams_twoLearnRates_alphaForget.bestParams(3)],'b','linewidth',2);
plot([1 2], [mdlCG15preS.fourParams_twoLearnRates_alphaForget.bestParams(3) mdlCG15effect.fourParams_twoLearnRates_alphaForget.bestParams(3)],'Color', [0.7 0 1],'linewidth',2);
xlim([0.75 2.25])
xticks([1 2])
xticklabels({'pre','post'})
title('alpha forget')

subplot(1,4,4); hold on;
plot([1 2], [mdlCG05pre.fourParams_twoLearnRates_alphaForget.bestParams(4) mdlCG05post.fourParams_twoLearnRates_alphaForget.bestParams(4)],'r','linewidth',2);
plot([1 2], [mdlCG14preS.fourParams_twoLearnRates_alphaForget.bestParams(4) mdlCG14effect.fourParams_twoLearnRates_alphaForget.bestParams(4)],'b','linewidth',2);
plot([1 2], [mdlCG15preS.fourParams_twoLearnRates_alphaForget.bestParams(4) mdlCG15effect.fourParams_twoLearnRates_alphaForget.bestParams(4)],'Color', [0.7 0 1],'linewidth',2);
xlim([0.75 2.25])
xticks([1 2])
xticklabels({'pre','post'})
title('beta')
legend('CG05', 'CG14', 'CG15')