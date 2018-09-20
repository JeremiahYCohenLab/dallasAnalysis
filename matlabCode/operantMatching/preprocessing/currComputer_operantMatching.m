function [root, sep] = currComputer_operantMatching()

if ismac
    root = '/Volumes/cooper/';
  %  root = '/Volumes/bbari1/';
    sep = '/';
elseif ispc
%     root = 'X:\';
    root = 'Y:\';
    sep = '\';
end