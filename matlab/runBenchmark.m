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
addpath(genpath(pwd));
PROBLEMS = {'BK1';'ex005';'Deb41';'Deb512a';'Deb512b';'Deb512c';'Deb513';
               'Deb521a';'Deb521b';'Deb53';'ZDT1';'ZDT2';'ZDT3';'ZDT4';
               'ZDT6';'DTLZ1';'DTLZ2';'DTLZ3';'DTLZ4';'DTLZ5';'DTLZ6';
               'DTLZ1n2';'DTLZ2n2';'DTLZ3n2';'DTLZ4n2';'DTLZ5n2';
               'DTLZ6n2';'Kursawe';'Fonseca';'L1ZDT4';'L2ZDT1';'L2ZDT2';
               'L2ZDT3';'L2ZDT4';'L2ZDT6';'L3ZDT1';'L3ZDT2';'L3ZDT3';
               'L3ZDT4';'L3ZDT6';'WFG1';'WFG2';'WFG3';'WFG4';'WFG5';
               'WFG6';'WFG7';'WFG8';'WFG9';'I1';'I2';'I3';'I4';'I5';
               'MOP1';'MOP2';'MOP3';'MOP4';'MOP5';'MOP6';'MOP7';
               'DPAM1';'DG01';'Far1';'FES1';'FES2';'FES3';'IKK1';'Jin1';
               'Jin2';'Jin3';'Jin4';'OKA1';'OKA2';'LRS1';'IM1';'LE1';
               'MHHM1';'MHHM2';'MLF1';'MLF2';'QV1';'Sch1';'SP1';'SSFYY1';
               'SSFYY2';'SK1';'SK2';'TKLY1';'VU1';'VU2';'VFM1';'ZLT1';
               'CL1';'lovison1';'lovison2';'lovison3';'lovison4';'lovison5';'lovison6'};
%PROBLEMS= {'Deb41'};
PROBLEMS_DIR = fullfile('..','problems');
ALGS_DIR = 'algs';         
EXP_DIR = fullfile('..','EXP_RESULTS');
TIM_DIR = fullfile('..','latex-template');
global nondomVectors;
global nondomTime;
global timeStamp;
if ~isdir(EXP_DIR)
    mkdir(EXP_DIR)
end

%% ======================================
% USER CONFIGURATIONS - CHANGE HERE:
EVAL_BUDGET_MULTIPLIER = 1000;% This value * DIM will be the evaluation budget for the corresponding problem
ALGS = {'MODIRECT'};
NUM_RUNS = [1] ;% for each alg in ALGS
%% ======================================

run_performance_experiments
run_timing_experiments
