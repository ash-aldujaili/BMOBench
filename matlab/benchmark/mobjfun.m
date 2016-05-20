function f=mobjfun(x, problemHandler, recordPareto)
% benchmarkFunc
%  evaluates the problem passed by problemHandler at x and keeps a record 
%  of the evaluation history by saving only new non-dominated points.
%   Takes:
%       x: an n-dimensional decision vector
%       problemHandler: handler for the multi-objective problem
%       recordPareto : a flag whether the Pareto front is recorded or not
%   Returns:
%       f: the m-dimensional objective vector f(x)
%--------------------------------------------------------------------------
%Copyright (c) 2016 by Abdullah Al-Dujaili
%
%This file is part of  the Multi Objective Competition & Benchmark
%It is free software: you can redistribute it and/or modify it under
%the terms of the GNU General Public License as published by the Free 
%Software Foundation, either version 3 of the License, or (at your option) 
%any later version.
%
%The associated files are distributed in the hope that it will be useful, but WITHOUT 
%ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
%FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
%more details.
%
%You should have received a copy of the GNU General Public License along.
%If not, see <http://www.gnu.org/licenses/>.

% Make sure x is a row-vector as problemHandler takes that
if (size(x,1) < size(x,2))
    x = x';
end

global nondomVectors;
global nondomTime;
global timeStamp;
% evaluate
timeStamp = timeStamp + 1;
f = problemHandler(x);
% decide to store or note
if (recordPareto)
    %pf = paretofront([ nondomVectors; f]);
    if isnondominated([ nondomVectors; f])%pf(end)
        nondomVectors(end+1,:) = f;
        nondomTime(end+1) = timeStamp;
    end
end

end