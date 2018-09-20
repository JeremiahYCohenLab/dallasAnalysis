function LH = likelihood(choice, probChoice)
% computes likelihood for learning models


Pf = choice.* probChoice;
nonzeros = Pf ~= 0;

if sum(sum(nonzeros)) < size(choice, 1) % if there are probabilities = 0
    LH = inf;
else
    Pf = Pf(nonzeros);
    LH = -1 * sum(log(Pf));
end