clc; clear; close all;

s = 'phasepack-matlab';

if isfolder(s)
   addpath(genpath(s))
else
    error('Please download the PhasePack package into this folder')
end

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2)

%% Success rate in terms of number of measurements m

genpar = struct;
genpar.d = 64;
genpar.signalType = 'gaussian'; % 'sinusoid'; % 

xitem = 'M';
xvalues = 2*genpar.d:genpar.d:7*genpar.d;  

%% Set up general parameters

params.verbose = true;            % verbose = wortreich -> more output details
params.runs = 20;                  % run several random trials for each scenario, and report average results
params.n = genpar.d;              % num of unknown elements
params.isComplex = false;         % use complex matrices? or just stick to real?
params.genpar = genpar;
params.successConstant = 0.05;

%% Benchmark algorithms

wf = struct('algorithm', 'wirtflow','initMethod', 'optimalspectral', 'label','Wirtinger Flow');
% fienup = struct('algorithm','fienup','initMethod', 'optimalspectral','label','Fienup');
gs = struct('algorithm','gerchbergsaxton', 'initMethod', 'optimalspectral', 'label','Gerchberg-Saxton');
pmax = struct( 'algorithm','phasemax', 'initMethod', 'optimalspectral', 'label','PhaseMax'); 
plift = struct('algorithm','phaselift', 'initMethod', 'optimalspectral', 'label','PhaseLift'); 
kv = struct('algorithm','custom', 'label','KV-CBO');          

algorithms = {wf, gs, pmax, plift, kv};

[finalResults, results] = runBenchmark(xitem, xvalues, algorithms, params);
