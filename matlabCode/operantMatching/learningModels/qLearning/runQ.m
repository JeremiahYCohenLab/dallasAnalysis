function [qFit, pe] = runQ(filename, testFlag)

[behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename);
behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);
[root, sep] = currComputer_operantMatching();

outcome = abs([behavStruct.allReward_R; behavStruct.allReward_L])';
choice = abs([behavStruct.allChoice_R; behavStruct.allChoice_L])';

% Fit all data
[qFit, pe] = qLearning_fitAll(choice, outcome, testFlag); sprintf('Q-learning finished')

BICs = [qFit.twoParams.BIC qFit.threeParams_alphaForget.BIC qFit.threeParams_bias.BIC qFit.fourParams.BIC];
if qFit.twoParams.BIC == min(BICs)
    pe = qFit.twoParams.pe;
elseif qFit.threeParams_alphaForget.BIC
    pe = qFit.threeParams_alphaForget.pe;
elseif qFit.threeParams_bias.BIC
    pe = qFit.threeParams_bias.pe;
elseif qFit.fourParams.BIC
    pe = qFit.fourParams.pe;
end
