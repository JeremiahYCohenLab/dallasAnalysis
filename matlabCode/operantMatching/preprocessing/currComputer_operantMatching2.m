function [root, sep] = currComputer_operantMatching2()

if ismac
    root = '/Users/Cooper/Desktop/';
  %  root = '/Volumes/bbari1/';
    sep = '/';
elseif ispc
     root = 'X:\';
%   root = 'Y:\';
    sep = '\';
end