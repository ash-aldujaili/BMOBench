global nondomVectors;
global nondomTime;
global timeStamp;
%% ======================================
% EXPERIMENT BODY : running ALGS on PROBLEMS NUM_RUNS times the budget
% a) Performance Evaluation
disp('Recording performance results...')
recordPareto = true;
for algIdx = 1 : numel(ALGS)
    alg = ALGS{algIdx};
    for problemIdx = 1: numel(PROBLEMS)
        problem = PROBLEMS{problemIdx};
        for run = 1 : NUM_RUNS(algIdx)
            disp(['    Solving ' problem ' by ' alg ', run:' num2str(run)])
            % initialize the problem:
            timeStamp = 0;
            nondomTime = [];
            nondomVectors = [];
            [~,l,u,~,~,~,~,m,~,~,~]=matampl(fullfile(PROBLEMS_DIR,[ problem '.nl']));
            dim = numel(l);
            %DIMENSION(end+1) = dim;
            numEvals = dim * EVAL_BUDGET_MULTIPLIER(algIdx);
            benchmarkFunc = @(x) mobjfun(x, @(x) arrayfun(@(y) matampl(x,y), 1 : m), recordPareto);
            %% ======================================
            % populate the algorithms of your choice here
            if strcmp(alg,'MODIRECT')
                MODIRECT(benchmarkFunc,l' , u', numEvals, m);			
            elseif strcmp(alg,'MORANDOM')
                MORANDOM(benchmarkFunc,l' , u', numEvals, m);
            else
                error('no such algorithm');
            end
            %% ======================================
            % record the output
            ofileName = sprintf('%s_%dD_%s_nfev%.1e_run%d.txt', problem, dim, alg, EVAL_BUDGET_MULTIPLIER(algIdx), run); 
            dlmwrite(fullfile(EXP_DIR,ofileName ), '# time stamp | objective vectors', 'delimiter',' ');
            dlmwrite(fullfile(EXP_DIR,ofileName ), [nondomTime(:), nondomVectors], 'delimiter',' ', '-append');       
        end
    end
end