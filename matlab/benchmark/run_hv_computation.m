for problemIdx = 1: numel(PROBLEMS)
    problem = PROBLEMS{problemIdx};
    disp([' Computing hv values for ', problem])
    % read roi
    roi = dlmread(fullfile(PROBLEMS_DIR,'roi',[problem '.bound']));
    roi_range = range(roi);
    normRefPoint = roi(2,:)./roi_range;
    % read problem specific
    [~,l,u,~,~,~,~,m,~,~,~]=matampl(fullfile(PROBLEMS_DIR,[ problem '.nl']));
    dim = numel(l);        
    for algIdx = 1 : numel(ALGS)
        alg = ALGS{algIdx};
        for run = 1 : NUM_RUNS(algIdx)
            %read pf data
            ifileName = sprintf('%s_%dD_%s_nfev%.1e_run%d.txt', problem, dim, alg, EVAL_BUDGET_MULTIPLIER(algIdx), run); 
            pf_data = csvread(fullfile(EXP_DIR,ifileName ),1,0);
            timeIdx = pf_data(:,1);
            normPf_data = bsxfun(@times,pf_data(:,2:end),1./roi_range);
            % Since computing hv is extremely intensive, we are sampling 20
            % points
            % uniformly within from the set of collected points
            % TODO find a really fast way to compute hv for FES3
            if strcmp(problem,'FES3') % 4-objective extremely computation expensive skip that
                timeIdx = timeIdx(1:2);
                incr_hv = zeros(numel(timeIdx),1);
                %disp('Skipping FES3 in this release because of computational complexity');
            else
                numSamples = 25;
                selectedTimeIdx = unique(int32 (linspace(1,numel(timeIdx),numSamples)));
                timeIdx = timeIdx(selectedTimeIdx);
                incr_hv = zeros(numel(timeIdx),1);
                % compute incremental hv
                
                for i = 1 : numel(timeIdx)
                    incr_hv(i) = hv(normPf_data(1:selectedTimeIdx(i),:)', normRefPoint);
                    %fprintf([num2str(i/numel(timeIdx)) '--']);
                end
                %fprintf('\n');
            end
            %% ======================================
            % record the output
            ofileName = sprintf('%s_%dD_%s_nfev%.1e_run%d_hv.txt', problem, dim, alg, EVAL_BUDGET_MULTIPLIER(algIdx), run); 
            dlmwrite(fullfile(EXP_DIR,ofileName ), '# time stamp | hv', 'delimiter','');
            dlmwrite(fullfile(EXP_DIR,ofileName ), [timeIdx(:), incr_hv], 'delimiter',' ', '-append');
        end
    end
end
