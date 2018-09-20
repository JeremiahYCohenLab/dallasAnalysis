function [preRwdHist, postRwdHist] = logRegByRwdHistComparison_opMD(animal)

%run function to generate lrm 
[preChange, preChoice,smoothWindows]= logRegByRwdHist_opMD('goodBehDays.xlsx',animal,'preS');  close;
[postChange, postChoice,~]= logRegByRwdHist_opMD('goodBehDays.xlsx',animal,'effect');  close;
tMax = length(smoothWindows);

%plot beta coeffs for single covariate type models
figure; suptitle(animal)
subplot(1,2,1); hold on
relevInds = 2:tMax+1;
coefVals = preChange.Coefficients.Estimate(relevInds);
CIbands = coefCI(preChange);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(smoothWindows,coefVals,errorL,errorU,'b','linewidth',2)

coefVals = postChange.Coefficients.Estimate(relevInds);
CIbands = coefCI(postChange);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(smoothWindows,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

title('Choice Switch on Reward History')
legend('pre', 'post')
xlabel('reward history window size')
ylabel('\beta Coefficient')
xlim([0 (smoothWindows(end) + 0.5)])


%plot beta coeffs for single covariate type models
subplot(1,2,2); hold on
relevInds = 2:tMax+1;
coefVals = preChoice.Coefficients.Estimate(relevInds);
CIbands = coefCI(preChoice);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(smoothWindows,coefVals,errorL,errorU,'b','linewidth',2)

coefVals = postChoice.Coefficients.Estimate(relevInds);
CIbands = coefCI(postChoice);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(smoothWindows,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

title('Choice on Side Reward History')
legend('pre', 'post')
xlabel('reward history window size')
ylabel('\beta Coefficient')
xlim([0 (smoothWindows(end) + 0.5)])