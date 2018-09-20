function model = actorCritic_fukai_fitAll(choice, outcome, testFlag)

if nargin < 3
    testFlag = 0;
end

% Initialize models
modelNames = {'base','base_bias','alphaLsame','alphaLsame_bias','full','full_bias'};
startValueCSVs = {'actorCritic_fukai_base_startValues.csv', ...
                  'actorCritic_fukai_base_bias_startValues.csv', ...
                  'actorCritic_fukai_alphaLsame_startValues.csv', ...
                  'actorCritic_fukai_alphaLsame_bias_startValues.csv', ...
                  'actorCritic_fukai_full_startValues.csv', ...
                  'actorCritic_fukai_full_bias_startValues.csv'};

% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alphaValue_range = [0 1];
alphaLearn_action_range = [0 1];
alphaForget_action_range = [0 1];
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
    if strcmp(modelNames{currMod}, 'base')
        b = [ alphaValue_range(2);   beta_range(2);  
             -alphaValue_range(1);  -beta_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@actorCritic_fukai_base, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actions] = ...
            actorCritic_fukai_base(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'base_bias')
        b = [ alphaValue_range(2);   beta_range(2);  bias_range(2);
             -alphaValue_range(1);  -beta_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@actorCritic_fukai_base_bias, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actions] = ...
            actorCritic_fukai_base_bias(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'alphaLsame')
        b = [ alphaValue_range(2);  alphaForget_action_range(2);  beta_range(2);
             -alphaValue_range(1); -alphaForget_action_range(1); -beta_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@actorCritic_fukai_alphaLsame, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actions] = ...
            actorCritic_fukai_alphaLsame(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'alphaLsame_bias')
        b = [ alphaValue_range(2);  alphaForget_action_range(2);  beta_range(2);  bias_range(2);
             -alphaValue_range(1); -alphaForget_action_range(1); -beta_range(1); -bias_range(1);];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@actorCritic_fukai_alphaLsame_bias, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actions] = ...
            actorCritic_fukai_alphaLsame_bias(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'full')
        b = [ alphaValue_range(2);  alphaLearn_action_range(2);  alphaForget_action_range(2);  beta_range(2); 
             -alphaValue_range(1); -alphaLearn_action_range(1); -alphaForget_action_range(1); -beta_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@actorCritic_fukai_full, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actions] = ...
            actorCritic_fukai_full(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'full_bias')
        b = [ alphaValue_range(2);  alphaLearn_action_range(2);  alphaForget_action_range(2);  beta_range(2);  bias_range(2);
             -alphaValue_range(1); -alphaLearn_action_range(1); -alphaForget_action_range(1); -beta_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@actorCritic_fukai_full_bias, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actions] = ...
            actorCritic_fukai_full_bias(model.(modelNames{currMod}).bestParams, choice, outcome);
    else
        error('Actor critic fukai: Model name not found')
    end
    
    model.(modelNames{currMod}).LH = -1 * LH(bestFit, :);
    model.(modelNames{currMod}).BIC = log(length(outcome))*numParam - 2*model.(modelNames{currMod}).LH;
    
    bestHess = squeeze(hess(bestFit, :, :));
    model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
    model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
end