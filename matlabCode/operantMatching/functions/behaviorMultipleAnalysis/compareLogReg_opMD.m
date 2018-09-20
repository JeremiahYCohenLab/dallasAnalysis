function [preAll, postAll] = compareLogReg_opMD(file, animal, pre, post, revForFlag, trialFlag)

if nargin < 6
    trialFlag = 0;
end
if nargin < 5
    revForFlag = 0;
end

%run function to generate lrm
if trialFlag
    [preAll, tMax]= combineLogReg_opMD(file, animal, pre, revForFlag);
    [postAll,~]= combineLogReg_opMD(file, animal, post, revForFlag);
else
    [preRwd, preNoRwd, preAll, s]= combineLogRegTime_opMD(file, animal, pre, revForFlag);
    [postRwd, postNoRwd, postAll,~]= combineLogRegTime_opMD(file, animal, post, revForFlag);
    tMax = s.tMax;
end


%plot beta coeffs for multiple covariate type model
figure;
subplot(1,2,1); hold on;
relevInds = 2:tMax+1;
coefVals = preAll.Coefficients.Estimate(relevInds);
CIbands = coefCI(preAll);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'b','linewidth',2)
else
    errorbar(((1:s.tMax)*s.binSize/1000),coefVals,errorL,errorU,'b','linewidth',2)
end

coefVals = postAll.Coefficients.Estimate(relevInds);
CIbands = coefCI(postAll);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlabel('Reward n trials back')
else
    errorbar(((1:s.tMax)*s.binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlim([0 (s.tMax*s.binSize/1000 + s.binSize/1000)])
    xlabel('Reward n seconds back')
end
title('Combined Model - Reward')
legend('pre', 'post')
ylabel('\beta Coefficient')


subplot(1,2,2); hold on;
relevInds = tMax+2:tMax*2+1;
coefVals = preAll.Coefficients.Estimate(relevInds);
CIbands = coefCI(preAll);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'b','linewidth',2)
else
    errorbar(((1:s.tMax)*s.binSize/1000),coefVals,errorL,errorU,'b','linewidth',2)
end

coefVals = postAll.Coefficients.Estimate(relevInds);
CIbands = coefCI(postAll);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlabel('Reward n trials back')
else
    errorbar(((1:s.tMax)*s.binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlim([0 (s.tMax*s.binSize/1000 + s.binSize/1000)])
    xlabel('Reward n seconds back')
end
title('Combined Model - Reward')
legend('pre', 'post')
ylabel('\beta Coefficient')

title('Combined Model - No Reward')
legend('pre', 'post')
ylabel('\beta Coefficient')


suptitle(animal)