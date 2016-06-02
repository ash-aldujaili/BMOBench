%--------------------------------------------------------------------------
% Compute the hv values of t




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
%--------------------------------------------------------------------------


PROBLEMS = {'CL1';'Deb512b';'Deb521b';'DTLZ1n2';'DTLZ3n2';'DTLZ5n2';'Far1';'Fonseca';'I4';'Jin1';'Jin3';'Kursawe';'L1ZDT4';'L3ZDT1';'L3ZDT2';'L3ZDT3';'L3ZDT4';'L3ZDT6';
  'VFM1';'VU1';'VU2';'ZDT1';'ZDT2';'ZDT3';'ZDT4';'ZDT6';'ZLT1';'Deb512a';'Deb521a';'Deb53';'DTLZ1';'DTLZ3';'DTLZ5';'ex005';'FES3';'I2';'I3';'IM1';'Jin4';'L2ZDT1';'L2ZDT2';'L2ZDT4';'L2ZDT6';'lovison1';
  'lovison3';'lovison5';'OKA1';'OKA2';'Sch1';'SK1';'SP1';'SSFYY2';'TKLY1';'WFG6';'WFG7';'WFG8';'BK1';'Deb41';'Deb512c';'DG01';'DTLZ2';'DTLZ4';'DTLZ6';'FES1';'I1';'I5';'L2ZDT3';'Jin2';'LE1';'lovison2';
  'lovison4';'lovison6';'MOP2';'MOP6';'QV1';'SK2';'SSFYY1';'WFG3';'MOP1';'MOP3';'MOP4';'MOP5';'MOP7';'MLF1';'MLF2';'Deb513';'DPAM1';'DTLZ2n2';'DTLZ4n2';'DTLZ6n2';'FES2';'LRS1';'MHHM1';'MHHM2';
  'WFG1';'WFG2';'WFG4';'WFG5';'WFG9';'IKK1'};
PROBLEMS_DIR = fullfile('..','problems');         
EXP_DIR = fullfile('..','EXP_RESULTS');


%% ======================================
% USER CONFIGURATIONS - CHANGE HERE:

ALGS = {'MORANDOM'};
EVAL_BUDGET_MULTIPLIER = [100];% This value * DIM will be the evaluation budget for the corresponding (problem,algorithm)
NUM_RUNS = [1] ;% for each alg in ALGS
%% ======================================

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
