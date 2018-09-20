figure;

subplot(1,3,1); hold on;
plot([1 2], [alphaNPE_CG05pre alphaNPE_CG05post],'r','linewidth',2);
plot([1 2], [alphaNPE_CG14preS alphaNPE_CG14effect],'b','linewidth',2);
plot([1 2], [alphaNPE_CG15preS alphaNPE_CG15effect],'Color', [0.7 0 1],'linewidth',2);
xlim([0.75 2.25])
title('alpha NPE')

subplot(1,3,2); hold on;
plot([1 2], [alphaPPE_CG05pre alphaPPE_CG05post],'r','linewidth',2);
plot([1 2], [alphaPPE_CG14preS alphaPPE_CG14effect],'b','linewidth',2);
plot([1 2], [alphaPPE_CG15preS alphaPPE_CG15effect],'Color', [0.7 0 1],'linewidth',2);
xlim([0.75 2.25])
title('alpha PPE')


subplot(1,3,3); hold on;
plot([1 2], [beta_CG05pre beta_CG05post], 'r','linewidth',2);
plot([1 2], [beta_CG14preS beta_CG14effect], 'b','linewidth',2);
plot([1 2], [beta_CG15preS beta_CG15effect],'Color', [0.7 0 1],'linewidth',2);
xlim([0.75 2.25])
title('beta')
legend('CG05', 'CG14', 'CG15')