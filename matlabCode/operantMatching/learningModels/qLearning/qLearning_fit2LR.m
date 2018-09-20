function [model] = qLearning_fit2LR(filename, testFlag, cellName)

if nargin < 2
    testFlag = 0;
    cellName = [];
end

[behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename);
behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);
[root, sep] = currComputer_operantMatching();

outcome = abs([behavStruct.allReward_R; behavStruct.allReward_L])';
choice = abs([behavStruct.allChoice_R; behavStruct.allChoice_L])';
ITI = [behavStruct.timeBtwn]';

if ~isempty(cellName)
    sessionName = regexp(filename, '.asc', 'split');
    spikeStruct = spikeProps_opMD(sessionName{1}, cellName);
    spikes = spikeStruct.maxCSrate';
    spikes = spikes/max(spikes);
end

% Initialize models
modelNames = {'fiveParams_twoLearnRates_alphaForget_bias'};
startValueCSVs = {'qLearningModel_5params_2learnRates_alphaForget_startValues.csv'};
              
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
alphaForget_s_range = [0 1];
alphaForget_l_range = [0 1];
tForget_range = [0 1];
beta_range = [0 double(intmax)];
v_range = [0 1];
ITIthresh_range = [0 120];
bias_range = [-5 5];

for currMod = 1:length(modelNames)
    if currMod == 4
        temp = load('C:\Users\cooper_PC\Desktop\GitHub repositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\qLearningModel_dual_startValues.mat');
        startValues = temp.startVals;
    else
        startValues = csvread(startValueCSVs{currMod});
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
    
    if strcmp(modelNames{currMod}, 'twoParams')
        b=[ alphaPPE_range(2);  beta_range(2);
           -alphaPPE_range(1); -beta_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_2params, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = qLearningModel_2params(model.(modelNames{currMod}).bestParams, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'fourParams_twoLearnRates_alphaForget')
        b=[ alphaNPE_range(2);  alphaPPE_range(2);  alphaForget_range(2); beta_range(2);
           -alphaNPE_range(1); -alphaPPE_range(1); -alphaForget_range(1); -beta_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_4params_2learnRates_alphaForget, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = ...
            qLearningModel_4params_2learnRates_alphaForget(model.(modelNames{currMod}).bestParams, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'fourParams_twoLearnRates_tForget')
        b=[ alphaNPE_range(2);  alphaPPE_range(2);  tForget_range(2); beta_range(2);
           -alphaNPE_range(1); -alphaPPE_range(1); -tForget_range(1); -beta_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_4params_2learnRates_tForget, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome, ITI);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = ...
            qLearningModel_4params_2learnRates_tForget(model.(modelNames{currMod}).bestParams, choice, outcome, ITI);
    end
    if strcmp(modelNames{currMod}, 'fourParams_twoLearnRates_alpha_spikes')
        b=[ alphaNPE_range(2);  alphaPPE_range(2);  tForget_range(2); beta_range(2);
           -alphaNPE_range(1); -alphaPPE_range(1); -tForget_range(1); -beta_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_4params_2learnRates_alpha_spikes, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome, ITI, spikes);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = ...
            qLearningModel_4params_2learnRates_alpha_spikes(model.(modelNames{currMod}).bestParams, choice, outcome, ITI, spikes);
    end
    if strcmp(modelNames{currMod}, 'fiveParams_twoLearnRates_alphaForget_bias')
        b=[ alphaNPE_range(2);  alphaPPE_range(2);  alphaForget_range(2); beta_range(2); bias_range(2);
           -alphaNPE_range(1); -alphaPPE_range(1); -alphaForget_range(1); -beta_range(1); -bias_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_5params_2learnRates_alphaForget, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = ...
            qLearningModel_5params_2learnRates_alphaForget(model.(modelNames{currMod}).bestParams, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'fiveParams_opponency')
        b=[ alphaNPE_range(2);  alphaPPE_range(2);  alphaForget_range(2); beta_range(2); v_range(2);
           -alphaNPE_range(1); -alphaPPE_range(1); -alphaForget_range(1); -beta_range(1); -v_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_5params_opponency, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe, model.(modelNames{currMod}).rBar] = ...
            qLearningModel_5params_opponency(model.(modelNames{currMod}).bestParams, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'fourParams_opponency_spikes')
        b=[ alphaNPE_range(2);  alphaPPE_range(2);  tForget_range(2); beta_range(2);
           -alphaNPE_range(1); -alphaPPE_range(1); -tForget_range(1); -beta_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_4params_opponency_spikes, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome, ITI, spikes);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = ...
            qLearningModel_4params_opponency_spikes(model.(modelNames{currMod}).bestParams, choice, outcome, ITI, spikes);
    end
    if strcmp(modelNames{currMod}, 'dual')
        b=[ alphaNPE_s_range(2); alphaPPE_s_range(2);  alphaForget_s_range(2); beta_range(2); ITIthresh_range(2);
            alphaNPE_l_range(2); alphaPPE_l_range(2);  alphaForget_l_range(2);
           -alphaNPE_s_range(1); -alphaPPE_s_range(1); -alphaForget_s_range(1); -beta_range(1); ITIthresh_range(1)
           -alphaNPE_l_range(1); -alphaPPE_l_range(1); -alphaForget_l_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_dual, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome, ITI);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q_s, model.(modelNames{currMod}).Q_l, model.(modelNames{currMod}).pe_s, model.(modelNames{currMod}).pe_l] = ...
            qLearningModel_dual(model.(modelNames{currMod}).bestParams, choice, outcome, ITI);
    end
    
    model.(modelNames{currMod}).LH = -1 * LH(bestFit, :);
    model.(modelNames{currMod}).BIC = log(length(outcome))*numParam - 2*model.(modelNames{currMod}).LH;
    
    bestHess = squeeze(hess(bestFit, :, :));
    model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
    model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
end
    
    
    