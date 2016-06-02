% This script verifies the C-coded problems against its AMPL-versions:
% this code is for development purposes only
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
%--------------------------------------------------------------------------

% SETUP - DO NOT CHANGE 
addpath(fullfile('benchmark'));
PROBLEMS = {'CL1';'Deb512b';'Deb521b';'DTLZ1n2';'DTLZ3n2';'DTLZ5n2';'Far1';'Fonseca';'I4';'Jin1';'Jin3';'Kursawe';'L1ZDT4';'L3ZDT1';'L3ZDT2';'L3ZDT3';'L3ZDT4';'L3ZDT6';
  'VFM1';'VU1';'VU2';'ZDT1';'ZDT2';'ZDT3';'ZDT4';'ZDT6';'ZLT1';'Deb512a';'Deb521a';'Deb53';'DTLZ1';'DTLZ3';'DTLZ5';'ex005';'FES3';'I2';'I3';'IM1';'Jin4';'L2ZDT1';'L2ZDT2';'L2ZDT4';'L2ZDT6';'lovison1';
  'lovison3';'lovison5';'OKA1';'OKA2';'Sch1';'SK1';'SP1';'SSFYY2';'TKLY1';'WFG6';'WFG7';'WFG8';'BK1';'Deb41';'Deb512c';'DG01';'DTLZ2';'DTLZ4';'DTLZ6';'FES1';'I1';'I5';'L2ZDT3';'Jin2';'LE1';'lovison2';
  'lovison4';'lovison6';'MOP2';'MOP6';'QV1';'SK2';'SSFYY1';'WFG3';'MOP1';'MOP3';'MOP4';'MOP5';'MOP7';'MLF1';'MLF2';'Deb513';'DPAM1';'DTLZ2n2';'DTLZ4n2';'DTLZ6n2';'FES2';'LRS1';'MHHM1';'MHHM2';
  'WFG1';'WFG2';'WFG4';'WFG5';'WFG9';'IKK1'};
PROBLEMS_DIR_AMPL = fullfile('..','problems');

tolsum = 0;
NUM_SAMPLES = 100;
% Loop over the problems and compare the values
for problemIdx = 1: numel(PROBLEMS)
        Y_VAL_AMPL = [];
        Y_VAL_C = [];
        problem = PROBLEMS{problemIdx};
        disp(['-Problem:' problem ])
        %% Ampl-coded settings
        [v, l1, u1, m1, y] = matc(problemIdx);
        [x,l,u,~,~,~,~,m,~,~,~]=matampl(fullfile(PROBLEMS_DIR_AMPL,[ problem '.nl']));
        if (sum(abs(l-l1'))>1e-5 || sum(abs(u-u1'))>1e-5)
            disp('LOWER');
            disp('AMPL');
            disp(l);
            disp('C');
            disp(l1);
            disp('UPPER');
            disp('AMPL');
            disp(u);
            disp('C');
            disp(u1);
            disp('mismatch in bounds/dimension');
        end
        dim = v;
        l = l1';
        u = u1';
        normalizedTestFunc =  @(x) arrayfun(@(y) matampl(x,y), 1 : m);
        testFunc = @(x) normalizedTestFunc((u-l).*x' + l);
        %% C-coded settings
        testFuncC = @(x) (mat2c(problemIdx,(u-l)'.*x + l'));
        %% Comparison code
        X_VAL = rand(NUM_SAMPLES, dim);
        for i = 1 : NUM_SAMPLES
            y = testFuncC(X_VAL(i,:));
            Y_VAL_AMPL = [Y_VAL_AMPL; testFunc(X_VAL(i,:))];
            Y_VAL_C = [Y_VAL_C; y]; 
        end
        % check for error
        tol = sum(sum(abs(Y_VAL_AMPL-Y_VAL_C)));
        tolsum = tolsum + tol;
        if tol> 1e-5 || isnan(tol)
            disp('Mismatch');
        end
end
disp(['Total absolute error ' num2str(tolsum)]);