% b) Timing Evaluation : this could have been incorporated in (a), but it
% would not be accurate as the nondominated vectors are being recorded and
% computed for every single evaluation budget which may lead in an increase
% in the algorithm's original runtime.
disp('Recording timing results...')
recordPareto = false;
timingMat = zeros(1,numel(ALGS));
feMat = zeros(1,numel(ALGS));
PROBLEMS = {'BK1','DPAM1','L3ZDT1','DTLZ3','FES3'};
for algIdx = 1 : numel(ALGS)
    alg = ALGS{algIdx};
    for problemIdx = 1: numel(PROBLEMS)
        problem = PROBLEMS{problemIdx};
        disp(['    Solving ' problem ' by ' alg])
        % initialize the problem:
        timeStamp = 0;
        [x,l,u,~,~,~,~,m,~,~,~]=matampl(fullfile(PROBLEMS_DIR,[ problem '.nl']));
        dim = numel(l);
        numEvals = dim * EVAL_BUDGET_MULTIPLIER;
        benchmarkFunc = @(x) mobjfun(x, @(x) arrayfun(@(y) matampl(x,y), 1 : m), recordPareto);
        %% ======================================
        tic;
        % populate the algorithms of your choice here
        if strcmp(alg,'MORANDOM')
            MORANDOM(benchmarkFunc,l' , u', numEvals, m);
        else
            error('no such algorithm');
        end
        %% ======================================
        timingMat(algIdx) = timingMat(algIdx) + toc;
        feMat(algIdx) = feMat(algIdx) + timeStamp;
    end
end
%% ======================================
% record timing data
timingMat = timingMat./feMat;
ofileName = sprintf('timing_%s.txt', strjoin(ALGS,'-'));
dlmwrite(fullfile(TIM_DIR,ofileName ), sprintf('%s',strjoin(ALGS,' ')), 'delimiter','');
dlmwrite(fullfile(TIM_DIR,ofileName ), timingMat, 'delimiter',' ', '-append'); 
%% ======================================

