function [model] = qLearning_fitAll(choice, outcome, testFlag)

if nargin < 3
    testFlag = 0;
end

% Initialize models
modelNames = {'twoParams','threeParams_alphaForget','threeParams_bias','fourParams'};
startValueCSVs = {'qLearningModel_2params_startValues.csv', 'qLearningModel_3params_alphaForget_startValues.csv', ...
                  'qLearningModel_3params_bias_startValues.csv', 'qLearningModel_4params_startValues.csv'};
              
% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alphaLearn_range = [0 1];
alphaForget_range = [0 1];
beta_range = [0 double(intmax)];
bias_range = [-5 5];

for currMod = 1:length(modelNames)
    startValues = csvread(startValueCSVs{currMod});
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
        b=[ alphaLearn_range(2);  beta_range(2);
           -alphaLearn_range(1); -beta_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_2params, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = qLearningModel_2params(model.(modelNames{currMod}).bestParams, choice, outcome);
        
    elseif strcmp(modelNames{currMod}, 'threeParams_alphaForget')
        b=[ alphaLearn_range(2);  alphaForget_range(2);  beta_range(2);
           -alphaLearn_range(1); -alphaForget_range(1); -beta_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_3params_alphaForget, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = qLearningModel_3params_alphaForget(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'threeParams_bias')
        b=[ alphaLearn_range(2);  beta_range(2);  bias_range(2);
           -alphaLearn_range(1); -beta_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_3params_bias, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = qLearningModel_3params_bias(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'fourParams')
        b=[ alphaLearn_range(2);  alphaForget_range(2);  beta_range(2);  bias_range(2);
           -alphaLearn_range(1); -alphaForget_range(1); -beta_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_4params, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = qLearningModel_4params(model.(modelNames{currMod}).bestParams, choice, outcome);
    else
        error('Q-learning: Model name not found')
    end
 
    model.(modelNames{currMod}).LH = -1 * LH(bestFit, :);
    model.(modelNames{currMod}).BIC = log(length(outcome))*numParam - 2*model.(modelNames{currMod}).LH;
    
    bestHess = squeeze(hess(bestFit, :, :));
    model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
    model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
%    model.(modelNames{currMod}).pe = pe(bestFit, :);
end