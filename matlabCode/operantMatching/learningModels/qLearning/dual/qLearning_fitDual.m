function [model] = qLearning_fitDual(filename, testFlag)

if nargin < 2
    testFlag = 0;
end

[behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename);
behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);
[root, sep] = currComputer_operantMatching();

outcome = abs([behavStruct.allReward_R; behavStruct.allReward_L])';
choice = abs([behavStruct.allChoice_R; behavStruct.allChoice_L])';
ITI = [behavStruct.timeBtwn]';


% Initialize models
modelNames = {'dual_tForget'};
startValueCSVs = {'qLearningModel_4params_2learnRates_alphaForget_startValues.csv'};
              
% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alphaNPE_range = [0 1];
alphaNPE_s_range = [0 1];
alphaNPE_l_range = [0 1];
alphaPPE_range = [0 1];
alphaPPE_s_range = [0 1];
alphaPPE_l_range = [0 1];
alphaForget_range = [0 1];
tForget_s_range = [0 100];
tForget_l_range = [0 100];
beta_range = [0 double(intmax)];
tau_range = [0.001 10];
Q_l_start_range = [0 1];

for currMod = 1:length(modelNames)
    if currMod == 1 
        temp = load('C:\Users\cooper_PC\Desktop\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\dual\qLearningModel_dual_tForget_startValues.mat');
        startValues = temp.startVals;
    elseif currMod == 999 
        temp = load('C:\Users\cooper_PC\Desktop\GitHub repositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\qLearningModel_dual_oneForget_startValues.mat');
        startValues = temp.startVals;
    elseif currMod == 999
        temp = load('C:\Users\cooper_PC\Desktop\GitHub repositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\qLearningModel_dual_tau_startValues.mat');
        startValues = temp.startVals;
    elseif currMod == 999
        temp = load('C:\Users\cooper_PC\Desktop\GitHub repositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\qLearningModel_dual_tau_noForget_startValues.mat');
        startValues = temp.startVals;
    elseif currMod == 999
        temp = load('C:\Users\cooper_PC\Desktop\GitHub repositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\qLearningModel_dual_oneForget_start_startValues.mat');
        startValues = temp.startVals;
    elseif currMod == 999
        temp = load('C:\Users\cooper_PC\Desktop\GitHub repositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\qLearningModel_dual_oneForget_2start_startValues.mat');
        startValues = temp.startVals;
    else
        temp = load('C:\Users\cooper_PC\Desktop\GitHub repositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\qLearningModel_fourParams_twoLearnRates_alphaForget_startValues.mat');
        startValues = temp.startVals;
    end
    if testFlag == 1
        startValues = startValues(1, :);
    end
   
    % initialize output variables
    runs = size(startValues, 1);
    allParams = zeros(size(startValues, 1), size(startValues, 2));
    LH = zeros(size(startValues, 1), 1);
    exitFl = zeros(size(startValues, 1), 1);
    hess = zeros(size(startValues, 1), size(startValues, 2), size(startValues, 2));
    numParam = size(startValues, 2);
    A=[eye(size(startValues, 2)); -eye(size(startValues, 2))];
    
    if strcmp(modelNames{currMod}, 'fourParams_twoLearnRates_alphaForget')
        b=[ alphaNPE_range(2);  alphaPPE_range(2);  alphaForget_range(2); beta_range(2);
           -alphaNPE_range(1); -alphaPPE_range(1); -alphaForget_range(1); -beta_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_4params_2learnRates_alphaForget, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = ...
            qLearningModel_4params_2learnRates_alphaForget(model.(modelNames{currMod}).bestParams, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'dual_tForget')
        b=[ alphaNPE_s_range(2); alphaPPE_s_range(2);  tForget_s_range(2);
            alphaNPE_l_range(2); alphaPPE_l_range(2);  tForget_l_range(2); beta_range(2);
           -alphaNPE_s_range(1); -alphaPPE_s_range(1); -tForget_s_range(1); 
           -alphaNPE_l_range(1); -alphaPPE_l_range(1); -tForget_l_range(1); -beta_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_dual_tForget, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome, ITI);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q_s, model.(modelNames{currMod}).Q_l, model.(modelNames{currMod}).pe_s, model.(modelNames{currMod}).pe_l] = ...
            qLearningModel_dual_tForget(model.(modelNames{currMod}).bestParams, choice, outcome, ITI);
    end
    if strcmp(modelNames{currMod}, 'dual_oneForget')
        b=[ alphaNPE_s_range(2); alphaPPE_s_range(2);  tForget_s_range(2);
            alphaNPE_l_range(2); alphaPPE_l_range(2);  beta_range(2);
           -alphaNPE_s_range(1); -alphaPPE_s_range(1); -tForget_s_range(1); 
           -alphaNPE_l_range(1); -alphaPPE_l_range(1); -beta_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_dual_oneForget, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome, ITI);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q_s, model.(modelNames{currMod}).Q_l, model.(modelNames{currMod}).pe_s, model.(modelNames{currMod}).pe_l] = ...
            qLearningModel_dual_oneForget(model.(modelNames{currMod}).bestParams, choice, outcome, ITI);
    end
    if strcmp(modelNames{currMod}, 'dual_oneForget_start')
        b=[ alphaNPE_s_range(2); alphaPPE_s_range(2);  tForget_s_range(2);
            alphaNPE_l_range(2); alphaPPE_l_range(2);  beta_range(2); Q_l_start_range(2);
           -alphaNPE_s_range(1); -alphaPPE_s_range(1); -tForget_s_range(1); 
           -alphaNPE_l_range(1); -alphaPPE_l_range(1); -beta_range(1); -Q_l_start_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_dual_oneForget_start, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome, ITI);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q_s, model.(modelNames{currMod}).Q_l, model.(modelNames{currMod}).pe_s, model.(modelNames{currMod}).pe_l] = ...
            qLearningModel_dual_oneForget_start(model.(modelNames{currMod}).bestParams, choice, outcome, ITI);
    end
    if strcmp(modelNames{currMod}, 'dual_oneForget_2start')
        b=[ alphaNPE_s_range(2); alphaPPE_s_range(2);  tForget_s_range(2);
            alphaNPE_l_range(2); alphaPPE_l_range(2);  beta_range(2); Q_l_start_range(2); Q_l_start_range(2);
           -alphaNPE_s_range(1); -alphaPPE_s_range(1); -tForget_s_range(1); 
           -alphaNPE_l_range(1); -alphaPPE_l_range(1); -beta_range(1); -Q_l_start_range(1); -Q_l_start_range(1)]; 
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_dual_oneForget_2start, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome, ITI);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q_s, model.(modelNames{currMod}).Q_l, model.(modelNames{currMod}).pe_s, model.(modelNames{currMod}).pe_l] = ...
            qLearningModel_dual_oneForget_2start(model.(modelNames{currMod}).bestParams, choice, outcome, ITI);
    end
    if strcmp(modelNames{currMod}, 'dual_tau')
        b=[ alphaNPE_s_range(2); alphaPPE_s_range(2);  tForget_s_range(2);
            alphaNPE_l_range(2); alphaPPE_l_range(2);  tForget_l_range(2); beta_range(2); tau_range(2);
           -alphaNPE_s_range(1); -alphaPPE_s_range(1); -tForget_s_range(1); 
           -alphaNPE_l_range(1); -alphaPPE_l_range(1); -tForget_l_range(1); -beta_range(1); -tau_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_dual_tau, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome, ITI);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q_s, model.(modelNames{currMod}).Q_l, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe_s, model.(modelNames{currMod}).pe_l] = ...
            qLearningModel_dual_tau(model.(modelNames{currMod}).bestParams, choice, outcome, ITI);
    end
    if strcmp(modelNames{currMod}, 'dual_tau_noForget')
        b=[ alphaNPE_s_range(2); alphaPPE_s_range(2); alphaNPE_l_range(2); alphaPPE_l_range(2);
            beta_range(2); tau_range(2);
           -alphaNPE_s_range(1); -alphaPPE_s_range(1); -alphaNPE_l_range(1); -alphaPPE_l_range(1);
           -beta_range(1); -tau_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_dual_tau_noForget, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome, ITI);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q_s, model.(modelNames{currMod}).Q_l, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe_s, model.(modelNames{currMod}).pe_l] = ...
            qLearningModel_dual_tau_noForget(model.(modelNames{currMod}).bestParams, choice, outcome, ITI);
    end
    
    model.(modelNames{currMod}).LH = -1 * LH(bestFit, :);
    model.(modelNames{currMod}).BIC = log(length(outcome))*numParam - 2*model.(modelNames{currMod}).LH;
    
    bestHess = squeeze(hess(bestFit, :, :));
    model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
    model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
end
    
    
    