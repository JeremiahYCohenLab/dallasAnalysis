%% R W D   H I S T   B E F O R E   S W I T C H
x = figure;
y = subplot(1,3,1)
hold on; 
histogram(rwdHistChangeCG14pre,'Normalization','probability','FaceColor',[0 0.7 1]); 
histogram(rwdHistChangeCG14post,'Normalization','probability','FaceColor', [0.7 0 1]);
xlabel('reward history before switch')
ylabel('probability')
legend('pre-lesion', 'post-lesion')
title('CG14')

yy = subplot(1,3,2)
hold on; 
histogram(rwdHistChangeCG15pre,'Normalization','probability','FaceColor',[0 0.7 1]); 
histogram(rwdHistChangeCG15post,'Normalization','probability','FaceColor', [0.7 0 1]);
xlabel('reward history before switch')
title('CG15')

yyy = subplot(1,3,3)
hold on; 
histogram(rwdHistChangeCG05pre,'Normalization','probability','FaceColor',[0 0.7 1]); 
histogram(rwdHistChangeCG05post,'Normalization','probability','FaceColor', [0.7 0 1]);
xlabel('reward history before switch')
title('CG05')

linkaxes([y yy yyy])
xlim([0 1])
set(x,'DefaultAxesFontName','Arial')


%% C O N S E C   N O   R W D S   B E F O R E   S W I T C H

x = figure;
y = subplot(1,3,1)
hold on; 
histogram(changeCG14pre,'Normalization','probability','FaceColor',[0 0.7 1]); 
histogram(changeCG14post,'Normalization','probability','FaceColor', [0.7 0 1]);
xlabel('consec no rewards before switch')
ylabel('probability')
legend('pre-lesion', 'post-lesion')
title('CG14')

yy = subplot(1,3,2)
hold on; 
histogram(changeCG15pre,'Normalization','probability','FaceColor',[0 0.7 1]); 
histogram(changeCG15post,'Normalization','probability','FaceColor', [0.7 0 1]);
xlabel('consec no rewards before switch')
title('CG15')

yyy = subplot(1,3,3)
hold on; 
histogram(changeCG05pre,'Normalization','probability','FaceColor',[0 0.7 1]); 
histogram(changeCG05post,'Normalization','probability','FaceColor', [0.7 0 1]);
xlabel('consec no rewards before switch')
title('CG05')

linkaxes([y yy yyy])
xlim([0 20])
set(x,'DefaultAxesFontName','Arial')

%% P R O B   S W I T C H   N O   R E W A R D
x = figure;
y = axes;
hold on;
plot([-40:length(behTblCG05.Prob_Switch)-41],behTblCG05.Prob_Switch, 'Color', [1 0 0], 'LineWidth',2)
plot([-33:length(behTblCG14.Prob_Switch)-34],behTblCG14.Prob_Switch, 'Color', [0 0.7 1], 'LineWidth',2)
plot([-38:length(behTblCG15.Prob_Switch)-39],behTblCG15.Prob_Switch, 'Color', [0.7 0 1], 'LineWidth',2)
vline(-0.5,'k')
xlabel('sessions')
ylabel('probability of switching after no reward')
legend('CG05','CG14','CG15')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])


%% P R O B   S T A Y   R E W A R D
x = figure;
y = axes;
hold on;
plot([-40:length(behTblCG05.Prob_Stay)-41],behTblCG05.Prob_Stay, 'Color', [1 0 0], 'LineWidth',2)
plot([-33:length(behTblCG14.Prob_Stay)-34],behTblCG14.Prob_Stay, 'Color', [0 0.7 1], 'LineWidth',2)
plot([-38:length(behTblCG15.Prob_Stay)-39],behTblCG15.Prob_Stay, 'Color', [0.7 0 1], 'LineWidth',2)
vline(-0.5,'k')
xlabel('sessions')
ylabel('probability of staying after reward')
legend('CG05','CG14','CG15')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])


%% F R A C T I O N    R E W A R D E D
x = figure;
y = axes;
hold on;
plot([-40:length(behTblCG05.Fraction_Rewarded)-41],behTblCG05.Fraction_Rewarded, 'Color', [1 0 0], 'LineWidth',2)
plot([-33:length(behTblCG14.Fraction_Rewarded)-34],behTblCG14.Fraction_Rewarded, 'Color', [0 0.7 1], 'LineWidth',2)
plot([-38:length(behTblCG15.Fraction_Rewarded)-39],behTblCG15.Fraction_Rewarded, 'Color', [0.7 0 1], 'LineWidth',2)
xlabel('sessions')
ylabel('fraction rewarded')
legend('CG05','CG14','CG15')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])
vline(-0.5,'k')
ylim([0.1 0.7])

%% F R A C T I O N    C O R R E C T
x = figure;
y = axes;
hold on;
plot([-40:length(behTblCG05.Fraction_Correct)-41],behTblCG05.Fraction_Correct, 'Color', [1 0 0], 'LineWidth',2)
plot([-33:length(behTblCG14.Fraction_Correct)-34],behTblCG14.Fraction_Correct, 'Color', [0 0.7 1], 'LineWidth',2)
plot([-38:length(behTblCG15.Fraction_Correct)-39],behTblCG15.Fraction_Correct, 'Color', [0.7 0 1], 'LineWidth',2)
vline(-0.5,'k')
xlabel('sessions')
ylabel('fraction correct')
legend('CG05','CG14','CG15')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])

%% N O R M   S W I T C H E S
x = figure;
y = axes;
hold on;
plot([-40:length(behTblCG05.Norm_Switches)-41],behTblCG05.Norm_Switches, 'Color', [1 0 0], 'LineWidth',2)
plot([-33:length(behTblCG14.Norm_Switches)-34],behTblCG14.Norm_Switches, 'Color', [0 0.7 1], 'LineWidth',2)
plot([-38:length(behTblCG15.Norm_Switches)-39],behTblCG15.Norm_Switches, 'Color', [0.7 0 1], 'LineWidth',2)
vline(-0.5,'k')
xlabel('sessions')
ylabel('normalized switches')
legend('CG05','CG14','CG15')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])


%% I T I   L I C K S
x = figure;
y = axes;
hold on;
plot([-40:length(behTblCG05.ITI_Licks)-41],behTblCG05.ITI_Licks, 'Color', [1 0 0], 'LineWidth',2)
plot([-33:length(behTblCG14.ITI_Licks)-34],behTblCG14.ITI_Licks, 'Color', [0 0.7 1], 'LineWidth',2)
plot([-38:length(behTblCG15.ITI_Licks)-39],behTblCG15.ITI_Licks, 'Color', [0.7 0 1], 'LineWidth',2)
vline(-0.5,'k')
xlabel('sessions')
ylabel('ITI licks / trial')
legend('CG05','CG14','CG15')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])

%% P R E - P O S T   S W I T C H

x = figure;
y = axes;
hold on;
plot([probSwitchCG05pre probSwitchCG05post], 'Color', [1 0 0], 'LineWidth',2)
plot([probSwitchCG14pre probSwitchCG14post], 'Color', [0 0.7 1], 'LineWidth',2)
plot([probSwitchCG15pre probSwitchCG15post], 'Color', [0.7 0 1], 'LineWidth',2)
ylabel('probability of switch given no reward')
legend('CG05','CG14','CG15')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])

%% P R E - P O S T   R W D   R A T E

x = figure;
y = axes;
hold on;
plot([rewardRateCG05pre rewardRateCG05post], 'Color', [1 0 0], 'LineWidth',2)
plot([rewardRateCG14pre rewardRateCG14post], 'Color', [0 0.7 1], 'LineWidth',2)
plot([rewardRateCG15pre rewardRateCG15post], 'Color', [0.7 0 1], 'LineWidth',2)
ylabel('reward rate')
legend('CG05','CG14','CG15')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])

%% P R E - P O S T   C O R R E C T   R A T E

x = figure;
y = axes;
hold on;
plot([correctRateCG05pre correctRateCG05post], 'Color', [1 0 0], 'LineWidth',2)
plot([correctRateCG14pre correctRateCG14post], 'Color', [0 0.7 1], 'LineWidth',2)
plot([correctRateCG15pre correctRateCG15post], 'Color', [0.7 0 1], 'LineWidth',2)
ylabel('correct rate')
legend('CG05','CG14','CG15')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])

%% P R E - P O S T   S T A Y

x = figure;
y = axes;
hold on;
plot([probStayCG05pre probStayCG05post], 'Color', [1 0 0], 'LineWidth',2)
plot([probStayCG14pre probStayCG14post], 'Color', [0 0.7 1], 'LineWidth',2)
plot([probStayCG15pre probStayCG15post], 'Color', [0.7 0 1], 'LineWidth',2)
ylabel('probability of Stay given no reward')
legend('CG05','CG14','CG15')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])

%% P R E - P O S T   A L P H A   N P E

x = figure;
y = axes;
hold on;
plot([alphaNPEcg05pre alphaNPEcg05post], 'Color', [1 0 0], 'LineWidth',2)
plot([alphaNPEcg14pre alphaNPEcg14post], 'Color', [0 0.7 1], 'LineWidth',2)
plot([alphaNPEcg15pre alphaNPEcg15post], 'Color', [0.7 0 1], 'LineWidth',2)
ylabel('NPE learning rate')
legend('CG05','CG14','CG15')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])

%% P R E - P O S T   A L P H A   P P E

x = figure;
y = axes;
hold on;
plot([alphaPPEcg05pre alphaPPEcg05post], 'Color', [1 0 0], 'LineWidth',2)
plot([alphaPPEcg14pre alphaPPEcg14post], 'Color', [0 0.7 1], 'LineWidth',2)
plot([alphaPPEcg15pre alphaPPEcg15post], 'Color', [0.7 0 1], 'LineWidth',2)
ylabel('PPE learning rate')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])

%% P R E - P O S T   B E T A

x = figure;
y = axes;
hold on;
plot([betaCG05pre betaCG05post], 'Color', [1 0 0], 'LineWidth',2)
plot([betaCG14pre betaCG14post], 'Color', [0 0.7 1], 'LineWidth',2)
plot([betaCG15pre betaCG15post], 'Color', [0.7 0 1], 'LineWidth',2)
ylabel('beta parameter')
set(x,'DefaultAxesFontName','Arial')
set(y, 'XTickLabel',[])
