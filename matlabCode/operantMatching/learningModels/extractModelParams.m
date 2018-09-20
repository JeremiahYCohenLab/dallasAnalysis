function extractModelParams

a = whos;
NPE = [];
PPE = [];
tForget = [];
beta = [];
v = [];
BIC = [];

for i = 1:length(a);
    NPE = [NPE eval(strcat(a(i).name, '.fiveParams_tForget_opponency.bestParams(1)'))];
    PPE = [PPE eval(strcat(a(i).name, '.fiveParams_tForget_opponency.bestParams(2)'))];
    tForget = [tForget eval(strcat(a(i).name, '.fiveParams_tForget_opponency.bestParams(3)'))];
    beta = [beta eval(strcat(a(i).name, '.fiveParams_tForget_opponency.bestParams(4)'))];
    v = [v eval(strcat(a(i).name, '.fiveParams_tForget_opponency.bestParams(5)'))];
    BIC = [BIC eval(strcat(a(i).name, '.fiveParams_tForget_opponency.BIC'))];
end
    