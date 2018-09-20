function [root, sep] = currComputer()

if ismac
    root = '/Volumes/cooper/';
  %  root = '/Volumes/bbari1/';
    sep = '/';
elseif ispc
%     root = 'X:\';
    root = 'Y:\';
    sep = '\';
end