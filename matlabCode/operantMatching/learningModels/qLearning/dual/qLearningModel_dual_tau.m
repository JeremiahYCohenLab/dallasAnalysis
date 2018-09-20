function [LH, probChoice, Q_s, Q_l, Q, pe_s, pe_l] = qLearningModel_dual_tau(startValues, choice, outcome, ITI)

alphaNPE_s = startValues(1);
alphaPPE_s = startValues(2);
tForget_s = startValues(3);
alphaNPE_l = startValues(4);
alphaPPE_l = startValues(5);
tForget_l = startValues(6);
beta = startValues(7);
tau = startValues(8);


trials = length(choice);
Q_s = zeros(trials,2);
Q_l = zeros(trials,2);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q_s(t+1, 2) = Q_s(t, 2);
        Q_l(t+1, 2) = Q_l(t, 2);
        pe_s(t) = outcome(t, 1) - Q_s(t, 1);
        pe_l(t) = outcome(t, 1) - Q_l(t, 1);
        if pe_s(t) < 0
            Q_s(t+1, 1) = Q_s(t, 1) + alphaNPE_s * pe_s(t);
        else
            Q_s(t+1, 1) = Q_s(t, 1) + alphaPPE_s * pe_s(t);
        end
        if pe_l(t) < 0
            Q_l(t+1, 1) = Q_l(t, 1) + alphaNPE_l * pe_l(t);
        else
            Q_l(t+1, 1) = Q_l(t, 1) + alphaPPE_l * pe_l(t);
        end
    else % left choice
        Q_s(t+1, 1) = Q_s(t, 1);
        Q_l(t+1, 1) = Q_l(t, 1);
        pe_s(t) = outcome(t, 2) - Q_s(t, 2);
        pe_l(t) = outcome(t, 2) - Q_l(t, 2);
        if pe_s(t) < 0
            Q_s(t+1, 2) = Q_s(t, 2) + alphaNPE_s * pe_s(t);
        else
            Q_s(t+1, 2) = Q_s(t, 2) + alphaPPE_s * pe_s(t);
        end
        if pe_l(t) < 0
            Q_l(t+1, 2) = Q_l(t, 2) + alphaNPE_l * pe_l(t);
        else
            Q_l(t+1, 2) = Q_l(t, 2) + alphaPPE_l * pe_l(t);
        end
    end
    Q_s(t+1, :) = Q_s(t+1, :)*exp(-tForget_s*ITI(t+1));
    Q_l(t+1, :) = Q_l(t+1, :)*exp(-tForget_l*ITI(t+1));
end

if choice(t, 1) == 1
    pe_s(trials) = outcome(end, 1) - Q_s(end, 1);
    pe_l(trials) = outcome(end, 1) - Q_l(end, 1);
else
    pe_s(trials) = outcome(end, 2) - Q_s(end, 2);
    pe_l(trials) = outcome(end, 2) - Q_l(end, 2);
end

% Find summed Q and call softmax  rule
Q_s = (1/tau)*ITI.*Q_s;
Q_l = tau*ITI.*Q_l;
Q = Q_s + Q_l;
probChoice = logistic([beta*(Q(:, 1)-Q(:, 2)), ...
                       beta*(Q(:, 2)-Q(:, 1))]);


% To calculate likelihood:
LH = likelihood_opMD(choice,probChoice);
end