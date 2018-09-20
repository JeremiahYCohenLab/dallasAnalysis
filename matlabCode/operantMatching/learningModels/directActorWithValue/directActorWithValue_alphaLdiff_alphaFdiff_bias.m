function [LH, probChoice, V, actionPreference] = directActorWithValue_alphaLdiff_alphaFdiff_bias(startValues, choice, outcome)

    
alphaLearn_action = startValues(1);
alphaLearn_value = startValues(2);
alphaForget_action = startValues(3);
alphaForget_value = startValues(4);
bias = startValues(5);

trials = length(choice);
actionPreference = zeros(trials,2);
probChoice = zeros(trials,2);
V = zeros(trials,1);


% Call learning rule
% initialize first probChoice
probChoice(1,:) = logistic([(actionPreference(1, 1) - actionPreference(1, 2)) + bias, ...
                            (actionPreference(1, 2) - actionPreference(1, 1)) - bias]);
for t = 1 : (trials-1)
    
    if choice(t, 1) == 1 % right choice
        V(t+1) = alphaForget_value*V(t) + alphaLearn_value*(outcome(t,1) - V(t)); % update global value
        actionPreference(t+1, 1) = alphaForget_action*actionPreference(t, 1) + alphaLearn_action*(1 - probChoice(t,1))*(outcome(t,1) - V(t+1));
        actionPreference(t+1, 2) = alphaForget_action*actionPreference(t, 2) - alphaLearn_action*probChoice(t,2)*(outcome(t,1) - V(t+1));
    else % left choice
       V(t+1) = alphaForget_value*V(t) + alphaLearn_value*(outcome(t,2) - V(t)); % update global value
        actionPreference(t+1, 2) = alphaForget_action*actionPreference(t, 2) + alphaLearn_action*(1 - probChoice(t,2))*(outcome(t,2) - V(t+1));
        actionPreference(t+1, 1) = alphaForget_action*actionPreference(t, 1) - alphaLearn_action*probChoice(t,1)*(outcome(t,2) - V(t+1));
    end
    probChoice(t+1, :) = logistic([(actionPreference(t+1,1) - actionPreference(t+1,2)) + bias, ...
                                   (actionPreference(t+1,2) - actionPreference(t+1,1)) - bias]);
end

% To calculate likelihood:
LH = likelihood(choice, probChoice);
end