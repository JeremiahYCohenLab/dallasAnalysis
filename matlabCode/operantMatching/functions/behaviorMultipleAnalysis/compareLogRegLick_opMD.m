function [preGLM, postGLM] = compareLogRegLick_opMD(file, animal, pre, post, revForFlag, trialFlag)

if nargin < 6
    trialFlag = 0;
end
if nargin < 5
    revForFlag = 0;
end

%run function to generate lrm 
if trialFlag
    [preGLM, preStayLickLat, preSwitchLickLat, tMax]= combineLogRegLickLat_opMD(file, animal, pre, 'revForFlag', revForFlag);
    [postGLM, postStayLickLat, postSwitchLickLat, ~]= combineLogRegLickLat_opMD(file, animal, post, 'revForFlag', revForFlag);
else
    [preGLM, preStayLickLat, preSwitchLickLat, binSize, timeMax]= combineLogRegLickLatTime_opMD(file, animal, pre, 'revForFlag', revForFlag);
    [postGLM, postStayLickLat, postSwitchLickLat, ~, ~]= combineLogRegLickLatTime_opMD(file, animal, post, 'revForFlag', revForFlag);
    timeBinEdges = [1000:binSize:timeMax];
    tMax = length(timeBinEdges) - 1;
end


%plot beta coeffs for multiple covariate type model
figure; hold on;
relevInds = 2:tMax+1;
coefVals = preGLM.Coefficients.Estimate(relevInds);
CIbands = coefCI(preGLM);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'b','linewidth',2)
else
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'b','linewidth',2)
end

coefVals = postGLM.Coefficients.Estimate(relevInds);
CIbands = coefCI(postGLM);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlabel('Reward n trials back')
else
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlim([0 (tMax*binSize/1000 + binSize/1000)])
    xlabel('Reward n seconds back')
end
title('LRM - Rewards on Licks')
suptitle(animal)
legend('pre', 'post')
ylabel('\beta Coefficient')


figure;
mag = [1 0 1];
cyan = [0 1 1];
set(gcf,'defaultAxesColorOrder',[mag; cyan]);

subplot(2,2,1)
yyaxis left; histogram(preStayLickLat, 30, 'Normalization', 'probability')
yyaxis right; histogram(postStayLickLat, 30, 'Normalization', 'probability')
legend('pre', 'post')
title('stay lick latency')

subplot(2,2,2)
yyaxis left; histogram(preSwitchLickLat, 30, 'Normalization', 'probability')
yyaxis right; histogram(postSwitchLickLat, 30, 'Normalization', 'probability')
legend('pre', 'post')
title('switch lick latency')

subplot(2,2,3)
yyaxis left; histogram(preStayLickLat, 30, 'Normalization', 'probability')
yyaxis right; histogram(preSwitchLickLat, 30, 'Normalization', 'probability')
legend('stay', 'switch')
title('pre lick latency')

subplot(2,2,4)
yyaxis left; histogram(postStayLickLat, 30, 'Normalization', 'probability')
yyaxis right; histogram(postSwitchLickLat, 30, 'Normalization', 'probability')
legend('stay', 'switch')
title('post lick latency')


suptitle(animal)