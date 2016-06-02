

fprintf('===Compiling mex files====\n');

fprintf('compiling paretofront.c...');
mex paretofront.c
fprintf('Done\n');

fprintf('compiling isnondominated.c...');
mex -I. isnondominated.c
fprintf('Done\n');

fprintf('compiling Hypervolume.cpp...');
mex -I. hv.cpp Hypervolume.cpp
fprintf('Done\n');

fprintf('compiling matc.c...');
mex -I. matc.c problems.c
fprintf('Done\n');

fprintf('===========================\n');

