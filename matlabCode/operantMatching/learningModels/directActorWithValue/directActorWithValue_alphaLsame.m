function [LH, probChoice, V, actionPreference] = directActorWithValue_alphaLsame(startValues, choice, outcome)

    
alphaLearn = startValues(1);

trials = length(choice);
actionPreference = zeros(trials,2);
probChoice = zeros(trials,2);
V = zeros(trials,1);


% Call learning rule
% initialize first probChoice
probChoice(1,:) = logistic([(actionPreference(1, 1) - actionPreference(1, 2)), ...
                            (actionPreference(1, 2) - actionPreference(1, 1))]);
for t = 1 : (trials-1)
    
    if choice(t, 1) == 1 % right choice
        V(t+1) = V(t) + alphaLearn*(outcome(t,1) - V(t)); % update global value
        actionPreference(t+1, 1) = actionPreference(t, 1) + alphaLearn*(1 - probChoice(t,1))*(outcome(t,1) - V(t+1));
        actionPreference(t+1, 2) = actionPreference(t, 2) - alphaLearn*probChoice(t,2)*(outcome(t,1) - V(t+1));
    else % left choice
       V(t+1) = V(t) + alphaLearn*(outcome(t,2) - V(t)); % update global value
        actionPreference(t+1, 2) = actionPreference(t, 2) + alphaLearn*(1 - probChoice(t,2))*(outcome(t,2) - V(t+1));
        actionPreference(t+1, 1) = actionPreference(t, 1) - alphaLearn*probChoice(t,1)*(outcome(t,2) - V(t+1));
    end
    probChoice(t+1, :) = logistic([(actionPreference(t+1,1) - actionPreference(t+1,2)), ...
                                   (actionPreference(t+1,2) - actionPreference(t+1,1))]);
end

% To calculate likelihood:
LH = likelihood(choice, probChoice);
end