function [LH, probChoice, V, actions] = actorCritic_fukai_alphaLsame_bias(startValues, choice, outcome)

alphaLearn = startValues(1);
alphaForget_action = startValues(2);
beta = startValues(3);
bias = startValues(4);

trials = length(choice);
V = zeros(trials,1);
actions = zeros(trials,2);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        V(t+1) = V(t) + alphaLearn*(outcome(t,1) - V(t)); % update global value
        
        actions(t+1, 1) = alphaForget_action*actions(t, 1) + alphaLearn * (outcome(t, 1) - V(t+1));
        actions(t+1, 2) = alphaForget_action*actions(t, 2);
        
    else % left choice
        V(t+1) = V(t) + alphaLearn*(outcome(t,2) - V(t));

        actions(t+1, 2) = alphaForget_action*actions(t, 2) + alphaLearn * (outcome(t, 2) - V(t+1));
        actions(t+1, 1) = alphaForget_action*actions(t, 1);

    end
end

% Call softmax rule
probChoice = logistic([beta*(actions(:,1) - actions(:,2)) + bias, ...
                       beta*(actions(:,2) - actions(:,1)) - bias]);

% To calculate likelihood:
LH = likelihood(choice, probChoice);
end