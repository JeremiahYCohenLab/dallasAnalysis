function model = dynamicThompsonSampling_fitAll(choice, outcome, testFlag, shortCircuitFlag)

if nargin < 3
    testFlag = 0;
    shortCircuitFlag = 0;
end

if shortCircuitFlag == 1 % break out of this function
    model.shortCircuit.bestParams = NaN;
    model.shortCircuit.LH = NaN;
    model.shortCircuit.BIC = NaN;
    model.shortCircuit.probChoice = NaN(length(choice), 2);
    return
end

% Initialize models
modelNames = {'sameAlphaLearn_noBias','sameAlphaLearn','full_noBias', 'full'};
startValueCSVs = {'dynamicThompsonSampling_sameAlphaLearn_noBias_startValues.csv', ...
                  'dynamicThompsonSampling_sameAlphaLearn_startValues.csv', ...
                  'dynamicThompsonSampling_full_noBias_startValues.csv', ...
                  'dynamicThompsonSampling_full_startValues.csv'};
              
% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alphaLearn_reward_range = [0 5];
alphaLearn_noreward_range = [0 5];
alphaForget_range = [0 1];
priorA_range = [0 100];
priorB_range = [0 100000];
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
    if strcmp(modelNames{currMod}, 'sameAlphaLearn_noBias')

        b=[ alphaLearn_reward_range(2);  alphaForget_range(2);  priorA_range(2);  priorB_range(2); 
           -alphaLearn_reward_range(1); -alphaForget_range(1); -priorA_range(1); -priorB_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@dynamicThompsonSampling_sameAlphaLearn_noBias, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).betaDist_right, model.(modelNames{currMod}).betaDist_left] = ...
            dynamicThompsonSampling_sameAlphaLearn_noBias(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'sameAlphaLearn')
        b=[ alphaLearn_reward_range(2);  alphaForget_range(2);  priorA_range(2);  priorB_range(2);  bias_range(2);
           -alphaLearn_reward_range(1); -alphaForget_range(1); -priorA_range(1); -priorB_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@dynamicThompsonSampling_sameAlphaLearn, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).betaDist_right, model.(modelNames{currMod}).betaDist_left] = ...
            dynamicThompsonSampling_sameAlphaLearn(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'full_noBias')
        b=[ alphaLearn_reward_range(2);  alphaLearn_noreward_range(2);  alphaForget_range(2);  priorA_range(2);  priorB_range(2); 
           -alphaLearn_reward_range(1); -alphaLearn_noreward_range(1); -alphaForget_range(1); -priorA_range(1); -priorB_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@dynamicThompsonSampling_full_noBias, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).betaDist_right, model.(modelNames{currMod}).betaDist_left] = ...
            dynamicThompsonSampling_full_noBias(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'full')
        b=[ alphaLearn_reward_range(2);  alphaLearn_noreward_range(2);  alphaForget_range(2);  priorA_range(2);  priorB_range(2);  bias_range(2); 
           -alphaLearn_reward_range(1); -alphaLearn_noreward_range(1); -alphaForget_range(1); -priorA_range(1); -priorB_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@dynamicThompsonSampling_full, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).betaDist_right, model.(modelNames{currMod}).betaDist_left] = ...
            dynamicThompsonSampling_full(model.(modelNames{currMod}).bestParams, choice, outcome);
    else 
        error('Dynamic Thompson sampling: Model name not found')
    end
    
    model.(modelNames{currMod}).LH = -1 * LH(bestFit, :);
    model.(modelNames{currMod}).BIC = log(length(outcome))*numParam - 2*model.(modelNames{currMod}).LH;
    
    bestHess = squeeze(hess(bestFit, :, :));
    model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
    model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
end