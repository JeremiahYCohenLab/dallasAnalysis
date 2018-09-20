function [LH, probChoice, X, A, S] = spModel_base(startValues, choice, outcome)

alphaLearn = startValues(1);
beta = startValues(2);

trials = length(choice);
X = zeros(trials, 2);
A = zeros(trials, 2);

consec_Rstack = 0;
consec_Lstack = 0;

X(1, :) = [1 1];

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice       
        A(t+1, 1) = A(t, 1) + alphaLearn*(outcome(t, 1)/X(t, 1) - A(t, 1));
        A(t+1, 2) = A(t, 2);
        
        consec_Lstack = consec_Lstack + 1;
        consec_Rstack = 0;
        
        X(t+1, 1) = 1;
        for i = 1:length(consec_Lstack)
            X(t+1, 2) = 1 + X(t, 2)*(1 - A(t, 2));
        end
    else % left choice
        A(t+1, 2) = A(t, 2) + alphaLearn*(outcome(t, 2)/X(t, 2) - A(t, 2));
        A(t+1, 1) = A(t, 1);
        
        consec_Rstack = consec_Rstack + 1;
        consec_Lstack = 0;
        
        X(t+1, 2) = 1;
        for i = 1:length(consec_Lstack)
            X(t+1, 1) = 1 + X(t, 1)*(1 - A(t, 1));
        end
    end
end
S = [A(:,1).*X(:,1) A(:,2).*X(:,2)];

% Call softmax  rule

probChoice = logistic([beta*(S(:, 1)-S(:, 2)), ...
                       beta*(S(:, 2)-S(:, 1))]);

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end