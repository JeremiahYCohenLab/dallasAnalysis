function model = directActorWithValue_fitAll(choice, outcome, testFlag)

if nargin < 3
    testFlag = 0;
end

% Initialize models
modelNames = {'alphaLsame', ...
              'alphaLsame_alphaFsame', ...
              'alphaLdiff_alphaFsame', ...
              'alphaLdiff_alphaFsame_bias', ...
              'alphaLdiff_alphaFdiff', ...
              'alphaLdiff_alphaFdiff_bias'};
startValueCSVs = {'directActorWithValue_alphaLsame.csv', ...
                  'directActorWithValue_alphaLsame_alphaFsame_startValues.csv', ...
                  'directActorWithValue_alphaLdiff_alphaFsame_startValues.csv', ...
                  'directActorWithValue_alphaLdiff_alphaFsame_bias_startValues.csv', ...
                  'directActorWithValue_alphaLdiff_alphaFdiff_startValues.csv', ...
                  'directActorWithValue_alphaLdiff_alphaFdiff_bias_startValues.csv'};
              
% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alphaLearn_action_range = [0 double(intmax)];
alphaLearn_value_range = [0 1];
alphaForget_action_range = [0 1];
alphaForget_value_range = [0 1];
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
    if strcmp(modelNames{currMod}, 'alphaLsame')
        b = [ alphaLearn_action_range(2);
             -alphaLearn_action_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@directActorWithValue_alphaLsame, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actionPreferences] = ...
            directActorWithValue_alphaLsame(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'alphaLsame_alphaFsame')
        b = [ alphaLearn_action_range(2);  alphaForget_action_range(2);
             -alphaLearn_action_range(1); -alphaForget_action_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@directActorWithValue_alphaLsame_alphaFsame, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actionPreferences] = ...
            directActorWithValue_alphaLsame_alphaFsame(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'alphaLdiff_alphaFsame')
        b = [ alphaLearn_action_range(2);  alphaLearn_value_range(2);  alphaForget_action_range(2);
             -alphaLearn_action_range(1); -alphaLearn_value_range(1); -alphaForget_action_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@directActorWithValue_alphaLdiff_alphaFsame, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actionPreferences] = ...
            directActorWithValue_alphaLdiff_alphaFsame(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'alphaLdiff_alphaFsame_bias')
        b = [ alphaLearn_action_range(2);  alphaLearn_value_range(2);  alphaForget_action_range(2);  bias_range(2);
             -alphaLearn_action_range(1); -alphaLearn_value_range(1); -alphaForget_action_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@directActorWithValue_alphaLdiff_alphaFsame_bias, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actionPreferences] = ...
            directActorWithValue_alphaLdiff_alphaFsame_bias(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'alphaLdiff_alphaFdiff')
        b = [ alphaLearn_action_range(2);  alphaLearn_value_range(2);  alphaForget_action_range(2);  alphaForget_value_range(2);
             -alphaLearn_action_range(1); -alphaLearn_value_range(1); -alphaForget_action_range(1); -alphaForget_value_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@directActorWithValue_alphaLdiff_alphaFdiff, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actionPreferences] = ...
            directActorWithValue_alphaLdiff_alphaFdiff(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'alphaLdiff_alphaFdiff_bias')
        b = [ alphaLearn_action_range(2);  alphaLearn_value_range(2);  alphaForget_action_range(2);  alphaForget_value_range(2);  bias_range(2);
             -alphaLearn_action_range(1); -alphaLearn_value_range(1); -alphaForget_action_range(1); -alphaForget_value_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@directActorWithValue_alphaLdiff_alphaFdiff_bias, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).V, model.(modelNames{currMod}).actionPreferences] = ...
            directActorWithValue_alphaLdiff_alphaFdiff_bias(model.(modelNames{currMod}).bestParams, choice, outcome);
    else 
        error('Direct actor with value: Model name not found')
    end
    
    model.(modelNames{currMod}).LH = -1 * LH(bestFit, :);
    model.(modelNames{currMod}).BIC = log(length(outcome))*numParam - 2*model.(modelNames{currMod}).LH;
    
    bestHess = squeeze(hess(bestFit, :, :));
    model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
    model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
end