function LH = likelihood_oMD(choice, probChoice)
% likelihood_oM    Computers likelihood for operant matching models
%   LH = likelihood_oM(choice, probChoice)
%   INPUTS
%       choice: trial x 2 matrix of choices
%           left column is left choices; right column is right choices
%           [1 0; 0 1; 0 1] corresponds to L->R->R
%       probChoice: trial x 2 matrix of probability of choices
%           left column is left choice probability; right column is right choice probability
%   OUTPUTS
%       LH: negative log likelihood

Pf = choice.* probChoice;
nonzeros = Pf ~= 0;

if sum(nonzeros(:)) < size(choice, 1) % if there are probabilities = 0, likely due to underflow
    LH = -1 * sum(log(0.01*ones(10*size(choice, 1), 1))); % give each choice a probability of 1% to give a bad likelihood
else
    Pf = Pf(nonzeros);
    LH = -1 * sum(log(Pf));
end