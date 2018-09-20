clear; close all; clc
gamma = linspace(0, 1, 1e3);
a = 1;
b = 1;
maxTrials = 1e3;
blockLength = [1e3];
rwdProb = [20 15];
p = DynamicForaging('RandomSeed',25,'BlockLength', blockLength,'MaxTrials', maxTrials, 'RewardProbabilities', rwdProb);
window = 1e3;
tF = 1e4;
% expKern = exp(-1/tF*(0:window));
expKern = ones(1, window+1);

rawData_plot = subplot(2, 1, 1); hold on
posterior_plot = subplot(2, 1, 2); hold on

posterior_Left =  gamma.^(a - 1).*(1 - gamma).^(b - 1);
posterior_Left = posterior_Left / sum(posterior_Left);
posterior_Right = gamma.^(a - 1).*(1 - gamma).^(b - 1);
posterior_Right = posterior_Right / sum(posterior_Right);

subplot(posterior_plot);
plot(posterior_Left)
plot(posterior_Right)
legend('Left', 'Right')
%%
iterA_L = a;
iterB_L = b;
iterA_R = a;
iterB_R = b; 
for currT = 1:p.MaxTrials
%     pause  

    Lsamp = betarnd(iterA_L, iterB_L);
    Rsamp = betarnd(iterA_R, iterB_R);
    
    Lsamp = iterA_L / (iterA_L + iterB_L);
    Rsamp = iterA_R / (iterA_R + iterB_R);
    
    if Lsamp > Rsamp
        p = p.inputChoice([1 0]);
    else
        p = p.inputChoice([0 1]);
    end
        
    if currT <= window
        minT = 1;
    else
        minT = currT - window;
    end
    
    rewardL = expKern(1:(currT - minT + 1))*p.AllRewards(currT:-1:minT, 1);
    rewardR = expKern(1:(currT - minT + 1))*p.AllRewards(currT:-1:minT, 2);
    choiceL = expKern(1:(currT - minT + 1))*p.AllChoices(currT:-1:minT, 1);
    choiceR = expKern(1:(currT - minT + 1))*p.AllChoices(currT:-1:minT, 2);
    
    posterior_Left =  gamma.^(rewardL + a - 1).*(1 - gamma).^(choiceL - rewardL + b - 1); 
    posterior_Left = posterior_Left / sum(posterior_Left);
    iterA_L = rewardL + a;
    iterB_L = choiceL - rewardL + b;
    posterior_Right = gamma.^(rewardR + a - 1).*(1 - gamma).^(choiceR - rewardR + b - 1);
    posterior_Right = posterior_Right / sum(posterior_Right);
    iterA_R = rewardR + a;
    iterB_R = choiceR - rewardR + b;
    
    subplot(rawData_plot); xlim([0 currT+1]); ylim([-1 1])
    if p.AllChoices(currT, 1) == 1 % left choice
        if p.AllRewards(currT, 1) == 1 % reward
            plot([currT, currT], [0 -1], 'k')
        else
            plot([currT, currT], [0 -0.5], 'k')
        end
    elseif p.AllChoices(currT, 2) == 1 % right choice
        if p.AllRewards(currT, 2) == 1 % reward
            plot([currT, currT], [0 1], 'k')
        else
            plot([currT, currT], [0 0.5], 'k')
        end
    end
    
    cla(posterior_plot)
    subplot(posterior_plot);
    plot(gamma, posterior_Left)
    plot(gamma, posterior_Right)
    legend('Left', 'Right')
end

%%
allC = p.AllChoices(:, 2);
allC(allC == 0) = -1;
allR = p.AllRewards(:, 2);
allR(p.AllRewards(:, 1) == 1) = -1;

tBack = 20;
rwdHx = NaN(length(allC), tBack);
choiceHx = NaN(length(allC), tBack);
for i = 1:tBack
    rwdHx(:, i) = [NaN(i, 1); allR(1:end - i)];
    choiceHx(:, i) = [NaN(i, 1); allC(1:end - i)];
end

rwdMod = fitglm([rwdHx ], allC == 1, 'distribution', 'binomial');
figure; 
subplot(3,1,1)
plot(rwdMod.Coefficients.Estimate(2:21))

%%
changeChoice = [false abs(diff(allC')) > 0];
changeHistogram = [];
changeHistogram_LtoR = [];
changeHistogram_RtoL = [];
for i = find(changeChoice == 1)
    if allC(i) == 1 % if a right lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allC(i-goBack) == -1 && allR(i-goBack) == 0 % if previous trial was a L lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram = [changeHistogram temp];
        changeHistogram_LtoR = [changeHistogram_LtoR temp];
    elseif allC(i) == -1 %if a left lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allC(i-goBack) == 1 && allR(i-goBack) == 0 % if previous trial was a R lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram = [changeHistogram temp];
        changeHistogram_RtoL = [changeHistogram_RtoL temp];
    end
end

subplot(3,1,2); hold on
histogram(changeHistogram_LtoR,0:max(changeHistogram)+2,'Normalization','probability')
histogram(changeHistogram_RtoL,0:max(changeHistogram)+2,'Normalization','probability')
legend('L -> R','R -> L')
xlabel('Consec No Rewards Before Switch (Trials)')

Rstays = [];
Lstays = [];
tempR = 0; tempL = 0;
for i = 2:length(allC)
    if allC(i) == 1 % right choice 
        tempR = tempR + 1;
    elseif allC(i) == -1 % left choice
        tempL = tempL + 1;
    end
    if allC(i) == 1 && allC(i-1) == -1 % right and previous was left
        Lstays = [Lstays tempL];
        tempL = 0;
    elseif allC(i) == -1 && allC(i-1) == 1 % left and previous was right
        Rstays = [Rstays tempR];
        tempR = 0;
    end 
end

subplot(3,1,3); hold on
histogram(Lstays,1:max([Rstays Lstays]),'Normalization','probability')
histogram(Rstays,1:max([Rstays Lstays]),'Normalization','probability')
xlabel('Stay Duration (Trials)')
legend('Left','Right')