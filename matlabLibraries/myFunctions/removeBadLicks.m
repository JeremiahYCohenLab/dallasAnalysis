function [licksOut] = removeBadLicks(licksIn, minimumLickTime)
% [licksOut] = removeBadLicks[licksIn, minimumLickTime]
% INPUT:
%   licksIn - vector of lick times
%   minimumLickTime - minimum time allowed between licks
% OUTPUT:
%   licksOut - vector of lick times

if nargin < 2
    minimumLickTime = 55;
end
removeInds = [0 diff([0 diff(licksIn)<minimumLickTime])] == 1;
while any(removeInds)
    licksIn(removeInds) = [];
    removeInds = [0 diff([0 diff(licksIn)<minimumLickTime])] == 1;
end
licksOut = licksIn;
end