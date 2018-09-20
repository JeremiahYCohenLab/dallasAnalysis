% clear; close all; clc
clear; clc

gamma = linspace(0, 1, 1e3);
a = 1;
b = 1;
alpha = 0.98;
maxTrials = 1e3;
blockLength = [50];
rwdProb = [40 10];
p = DynamicForaging('RandomSeed',25,'BlockLength', blockLength,'MaxTrials', maxTrials, 'RewardProbabilities', rwdProb);
p0 = gamma.^(a - 1).*(1 - gamma).^(b - 1); p0 = p0/sum(p0);
iterposterior_Left =  p0; iterposterior_Left = iterposterior_Left/sum(iterposterior_Left);
iterposterior_Right = p0; iterposterior_Right = iterposterior_Right/sum(iterposterior_Right);
% posterior_Left = p0; posterior_Left = posterior_Left/sum(posterior_Left);
% posterior_Right = p0; posterior_Right = posterior_Right/sum(posterior_Right);

figure
rawData_plot = subplot(2, 1, 1); hold on
iterposterior_plot = subplot(2, 1, 2); hold on; title('Iterative Posterior / Iterative Prior')
plot(gamma, iterposterior_Left)
plot(gamma, iterposterior_Right)
legend('Left', 'Right')

for currT = 1:p.MaxTrials
    pause
    
%     % Thompson sample
%     scaledLeft = round(iterposterior_Left*1e6);
%     rnd_Left = round(rand*sum(round(iterposterior_Left*1e6)));
%     cumsum_scaledLeft = cumsum(scaledLeft); 
%     tmp = cumsum_scaledLeft - rnd_Left;
%     tmp(tmp < 0) = -1;
%     tmp(tmp > 0) = 1;
%     myInd = find([0 diff(tmp)] ~= 0); 
%     if length(myInd) > 1
%         myInd = myInd(1);
%     elseif isempty(myInd)
%         myInd = 1;
%     end
%     sampleLeft = gamma(myInd);
%     
%     scaledRight = round(iterposterior_Right*1e6);
%     rnd_Right = round(rand*sum(round(iterposterior_Right*1e6)));
%     cumsum_scaledRight = cumsum(scaledRight);
%     tmp = cumsum_scaledRight - rnd_Right;
%     tmp(tmp < 0) = -1;
%     tmp(tmp > 0) = 1;
%     myInd = find([0 diff(tmp)] ~= 0); 
%     if length(myInd) > 1
%         myInd = myInd(1);
%     elseif isempty(myInd)
%         myInd = 1;
%     end
%     sampleRight = gamma(myInd);

%     % Select highest mean
    sampleLeft = sum(gamma.*iterposterior_Left);
    sampleRight = sum(gamma.*iterposterior_Right);
    
    if rand < 0.05
        if sampleLeft > sampleRight
            sampleRight = sampleRight + 1;
        else
            sampleLeft = sampleLeft + 1;
        end
    end
    
    
    % make a choice
    if sampleLeft > sampleRight
        p = p.inputChoice([1 0]);
        r = sum(p.AllRewards(p.Trial, :));
        likelihood = (gamma.^r).*(1 - gamma).^(1 - r);

        iterposterior_Left = likelihood.*iterposterior_Left;
        iterposterior_Left = iterposterior_Left/sum(iterposterior_Left);
    else
        p = p.inputChoice([0 1]);
        r = sum(p.AllRewards(p.Trial, :));
        likelihood = (gamma.^r).*(1 - gamma).^(1 - r);
        
        iterposterior_Right = likelihood.*iterposterior_Right;
        iterposterior_Right = iterposterior_Right/sum(iterposterior_Right);
    end
    
    iterposterior_Left = alpha*iterposterior_Left + (1 - alpha)*p0;
    iterposterior_Left = iterposterior_Left/sum(iterposterior_Left);
    
    iterposterior_Right = alpha*iterposterior_Right + (1 - alpha)*p0;
    iterposterior_Right = iterposterior_Right/sum(iterposterior_Right);

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
    
    cla(iterposterior_plot)
    subplot(iterposterior_plot)
    plot(gamma, iterposterior_Left)
    plot(gamma, iterposterior_Right)
    legend('Left','Right')
    
    P_XgY = 0;
    for i = 1:length(gamma)
        P_XgY = P_XgY + sum(iterposterior_Left(1:i))*iterposterior_Right(i);
    end
    title(sprintf('Probability of right choice %d%%', round(P_XgY*100)))
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

rwdMod = fitglm([rwdHx], allC == 1, 'distribution', 'binomial');
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