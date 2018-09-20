%% VI, VI
failures = 0:100;
pr = 0.4;
pl = 0.1;


prt = 1 - (1 - pr).^(failures+1);
plt = 1 - (1 - pl).^(failures+1);
% prt = prt / sum(prt);

allPR = 0:.01:1;
r_R = NaN(size(allPR));
r_L = NaN(size(allPR));
% allP1 = 0.9;
for i = 1:length(allPR)
    pt = (1 - allPR(i)).^(failures).*allPR(i);
    pt = pt / sum(pt);
    
    r_R(i) = sum(prt.*pt);
    r_L(end+1-i) = sum(plt.*pt);
end

% r_L = pl * ones(size(r_R));
r_full = r_R.*allPR + r_L.*(1-allPR);

figure; hold on;
plot(allPR,r_R)
plot(allPR,r_L)
plot(allPR,r_full);
[~,maxInd] = max(r_full);
plot(allPR(maxInd),r_full(maxInd),'ro')

legend('<r|R>','<r|L>','<r>','max <r>')
title(sprintf('Right: %s\nLeft: %s',num2str(pr),num2str(pl)))

xlabel('P(R)')
ylabel('Reward')

%% VI (right), VR (left)
failures = 0:100;
pr = 0.3;  % VI
pl = 0.35; % VR


prt = 1 - (1 - pr).^(failures+1);
plt = pl; %1 - (1 - pl).^(failures+1);
% prt = prt / sum(prt);

allPR = 0:.01:1;
r_R = NaN(size(allPR));
r_L = NaN(size(allPR));
% allP1 = 0.9;
for i = 1:length(allPR)
    pt = (1 - allPR(i)).^(failures).*allPR(i);
    pt = pt / sum(pt);
    
    r_R(i) = sum(prt.*pt);
    r_L(end+1-i) = sum(plt.*pt);
end

r_full = r_R.*allPR + r_L.*(1-allPR);

figure; hold on;
plot(allPR,r_R)
plot(allPR,r_L)
plot(allPR,r_full);
[~,maxInd] = max(r_full);
plot(allPR(maxInd),r_full(maxInd),'ro')

legend('<r|R>','<r|L>','<r>','max <r>')
title(sprintf('Right: %s\nLeft: %s',num2str(pr),num2str(pl)))

xlabel('P(R)')
ylabel('Reward')

%% VI (right), VR (left)
failures = 0:100;
pl = 0.4;
pr = 0:.01:pl;

matchMaxDistance = NaN(size(pr));

for currR = 1:length(pr)
    prt = 1 - (1 - pr(currR)).^(failures+1);
    plt = pl;

    allPR = 0:.001:1;
    r_R = NaN(size(allPR));
    r_L = NaN(size(allPR));
    % allP1 = 0.9;
    for i = 1:length(allPR)
        pt = (1 - allPR(i)).^(failures).*allPR(i);
        pt = pt / sum(pt);

        r_R(i) = sum(prt.*pt);
        r_L(end+1-i) = sum(plt.*pt);
    end

    r_full = r_R.*allPR + r_L.*(1-allPR);

    [~,maxInd] = max(r_full);
    [~, matchInd] = min(abs(r_L - r_R));
    matchMaxDistance(currR) = (allPR(matchInd) - allPR(maxInd));
end

figure; plot(pr,matchMaxDistance)
xlabel('Right reward (baited)')
ylabel('Distance in %')
title(sprintf('Left: %s',num2str(pl)))
%% Arbitrary task
i1 = .8; s1 = -.6; i2 = 0.4; s2 = 0;
match = (i2-i1)/(s1-s2)
maximize = (i2-i1-s2)/(2*s1-2*s2)

%% Competitive Foraging Task
a_R = 0.35; % add at this percent
u_R = 0.7; % remove at this percent
a_L = 0.25;
u_L = 0.1;

failures = 1:1000;
pr_R = (a_R * (1 - (1 - a_R).^failures .* (1 - u_R).^failures)) ./ (1 - (1 - a_R) .* (1 - u_R)); % P(reward|t trails since choosing)
pr_L = (a_L * (1 - (1 - a_L).^failures .* (1 - u_L).^failures)) ./ (1 - (1 - a_L) .* (1 - u_L));

probabilitySet = 0:.01:1;
r_R = NaN(size(probabilitySet)); 
r_L = NaN(size(probabilitySet));
for i = 1:length(probabilitySet)
    
    pt = (1 - probabilitySet(i)).^(failures).*probabilitySet(i); % P(t trials since choosing)
    pt = pt / sum(pt); % normalize
    
    r_R(i) = sum(pr_R.*pt); % mean reward on R side
    r_L(end+1-i) = sum(pr_L.*pt); % mean reward on L side
end

r_full = r_R.*probabilitySet + r_L.*(1-probabilitySet); % mean reward total

figure; hold on;
plot(probabilitySet,r_R)
plot(probabilitySet,r_L)
plot(probabilitySet,r_full);
[~,maxInd] = max(r_full);
plot(probabilitySet(maxInd),r_full(maxInd),'ro') % overall maximum

legend('<r|R>','<r|L>','<r>','max <r>')
title(sprintf('Right: lambda - %s, mu - %s\nLeft: lambda - %s, mu -%s',num2str(a_R),num2str(u_R),num2str(a_L),num2str(u_L)))
xlabel('P(R)')
ylabel('Reward')