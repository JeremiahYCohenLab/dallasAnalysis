function [glm_rwd, glm_rwdTime] = compareLogRegModels_opMD(xlFile, animal, category)


[~, combinedRewardsMatx, glm_rwd, trials] = combineLogReg_opMD(xlFile, animal, category);      %for trial-wise
[aicTrial, bicTrial] = aicbic(glm_rwd.LogLikelihood, trials, sum(~isnan(combinedRewardsMatx(1,:))));

[glm_rwdTime, ~, ~, binSize, bins] = logRegTime_opMD(xlFile, animal, category);      %for time-wise
[aicTime, bicTime] = aicbic(glm_rwdTime.LogLikelihood, bins, sum(~isnan(combinedRewardsMatx(1,:))));

figure;
subplot(1,3,1)
plot([1 2], [aicTrial aicTime], 'b') 
ylabel('AIC score')
xlim([0 3])
xticks([1 2])
xticklabels({'trial', 'time'})

subplot(1,3,2) 
plot([1 2], [bicTrial bicTime], 'b') 
ylabel('BIC score')
xlim([0 3])
xticks([1 2])
xticklabels({'trial', 'time'})


subplot(1,3,3); hold on;
relevInds = 2:21;
coefVals = glm_rwdTime.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwdTime);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar([1:20],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

relevInds = 2:21;
coefVals = glm_rwd.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwd);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar([1:20],coefVals,errorL,errorU,'b','linewidth',2)

timeLeg = [num2str(binSize/1000) 's time bins']; 
legend(timeLeg, 'trials')
xlabel('Outcome n trials/bins Back')
ylabel('\beta Coefficient')
xlim([0 21])
suptitle([animal ' ' category])