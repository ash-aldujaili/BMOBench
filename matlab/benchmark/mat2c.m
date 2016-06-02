function [ y] = mat2c(problem, x)
%MAT2C Wrapper for matc coded problems
% INPUT problem : problem's name as a string
%       x       : 1xn variable
[~, ~, ~, ~, y]=matc(problem,x);
end

