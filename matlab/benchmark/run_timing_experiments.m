% b) Timing Evaluation : this could have been incorporated in (a), but it
% would not be accurate as the nondominated vectors are being recorded and
% computed for every single evaluation budget which may lead in an increase
% in the algorithm's original runtime.
global timeStamp;
disp('Recording timing results...')
BUDGET_MULTIPLIER = [10000,1000,100, 10];
recordPareto = false;
timingMat = zeros(numel(BUDGET_MULTIPLIER),numel(ALGS));
feMat = zeros(numel(BUDGET_MULTIPLIER),numel(ALGS));
PROBLEMS = {'BK1','DPAM1','L3ZDT1','DTLZ3','FES3'};

for budMultIdx = 1: numel(BUDGET_MULTIPLIER)
    for algIdx = 1 : numel(ALGS)
        alg = ALGS{algIdx};
        for problemIdx = 1: numel(PROBLEMS)
            problem = PROBLEMS{problemIdx};
            disp(['    Solving ' problem ' by ' alg])
            % initialize the problem:
            timeStamp = 0;
            [~,l,u,~,~,~,~,m,~,~,~]=matampl(fullfile(PROBLEMS_DIR,[ problem '.nl']));
            dim = numel(l);
            numEvals = BUDGET_MULTIPLIER(budMultIdx);
            benchmarkFunc = @(x) mobjfun(x, @(x) arrayfun(@(y) matampl(x,y), 1 : m), recordPareto);
            %% ======================================
            tic;
            % populate the algorithms of your choice here
            if strcmp(alg,'MODIRECT')
                MODIRECT(benchmarkFunc,l' , u', numEvals, m);
            elseif strcmp(alg,'MORANDOM')
                MORANDOM(benchmarkFunc,l' , u', numEvals, m);
            else
                error('no such algorithm');
            end
            %% ======================================
            timingMat(budMultIdx,algIdx) = timingMat(budMultIdx,algIdx) + toc;
            feMat(budMultIdx,algIdx) = feMat(budMultIdx,algIdx) + timeStamp;
        end
    end
end
%% ======================================
% record timing data
timingMat = timingMat./feMat;
ofileName = sprintf('timing.txt');
dlmwrite(fullfile(TIM_DIR,ofileName ), sprintf('%s',strjoin(['x' ALGS],' ')), 'delimiter','');
dlmwrite(fullfile(TIM_DIR,ofileName ), [BUDGET_MULTIPLIER',timingMat], 'delimiter',' ', '-append'); 
%% ======================================

