%{
Author: Amine Laghaout
Date: 2015-01-30

Script for a batch run of Reporter() over several workers.

====
delete(gcp('nocreate')); clear; clc; close all; delete *.mat; 
myJob = batch('BatchRun', 'Profile', 'local', 'Pool', 16);  % G-bar local
myJob = batch('BatchRun', 'Profile', 'DTU_Torque', 'Pool', 31);  % DTU_Torque
myJob = batch('BatchRun', 'Profile', 'local', 'Pool', 1);   % Laptop
====
%}
function BatchRun()

% ALAWAYS CHECK THAT...
% - `iter' = 5 in AdaptiveDecision()
% - `S' = myPool.numWorkers in main()
% - `S' = 1 if we're calling mainGlobal() because SEQ not NaN
% - `SEQ' = NaN in Reporter to execute adaptive feedback
% - Number of workers on DTU_Torque should be set to 32
% - `ImTol' = 1e-6 in Distinguishability()

%% Produce figures from reporters
%{
matfile = cellstr(['20141222T151417'; '20141224T204741'; '20141223T185621'; '20141225T024522'; '20141223T201836'; '20141227T140928'; '20141224T180934'; '20141224T204426'; '20150205T140036']);

for j = 1:length(matfile)
    PrintReport(['Reporter', matfile{j}]);
end
%}

%% Compare local optimzations with those recorded globally
%{
matfile = cellstr(['20141227T174438'; '20141228T164245'; '20141229T021302'; '20141229T210728'; '20141229T210953'; '20141229T212420'; '20141229T212539'; '20141229T213314'; '20141229T215323'; '20141230T225734'; '20141230T231105']);
FoM = cellstr(['B         '; 'D         '; 'DiscErr   '; 'avgOverlap'; 'Bpi       ']);

for j = 1:length(matfile)
    
    matfilePath = ['./reports/GlobalReporter', matfile{j}, '.mat'];
    if exist(matfilePath, 'file') == 2
        load(matfilePath);
    else
        error('%s.mat does not exist.', matfile)
    end

    SEQ = [B.params; D.params; DiscErr.params; avgOverlap.params];
    fprintf('%s:\n', matfile{j});
    fprintf('F%d M%d N%d C%d %s %s [%.2f %.2f %d] eta=%.2f\n', FockDim, M, N, C, func2str(Operation), func2str(Pi), ParamsBounds(1), ParamsBounds(2), ParamsBounds(3), eta);
    fprintf('---global:\n')
    for i = 1:4
        [D, B, ~, ~, avgOverlap, DiscErr, ~] = ...
            main(FockDim, Operation, ParamsBounds, Pi, eta, M, C, N, SEQ(i, :), 'D');
        fprintf('%.4f\t%.4f\t%.4f\t%.4f\t[', B, D, DiscErr, avgOverlap);
        fprintf('%.3f ', SEQ(i, :));
        fprintf(']\n');
    end
    
    for k = 1:length(FoM)
        fprintf('---local over %s:\n', FoM{k});
        [D, B, ~, ~, avgOverlap, DiscErr, ~] = ...
            main(FockDim, Operation, ParamsBounds, Pi, eta, M, C, N, NaN, FoM{k});
        fprintf('%.4f\t%.4f\t%.4f\t%.4f\n', B, D, DiscErr, avgOverlap);    
    end
    
    fprintf('\n');
    
end
%}

%% Compare global optimzations
%{
matfile = cellstr(['20141227T174438'; '20141228T164245'; '20141229T021302'; '20141229T210728'; '20141229T210953'; '20141229T212420'; '20141229T212539'; '20141229T213314'; '20141229T215323'; '20141230T225734'; '20141230T231105']);
FoM = 'D';

for j = 1:length(matfile)
    
    matfilePath = ['./reports/GlobalReporter', matfile{j}, '.mat'];
    if exist(matfilePath, 'file') == 2
        load(matfilePath);
    else
        error('%s.mat does not exist.', matfile)
    end

    SEQ = [B.params; D.params; DiscErr.params; avgOverlap.params];
    fprintf('%s:\n', matfile{j});
    fprintf('F%d M%d N%d C%d %s %s [%.2f %.2f %d] eta=%.2f\n', FockDim, M, N, C, func2str(Operation), func2str(Pi), ParamsBounds(1), ParamsBounds(2), ParamsBounds(3), eta);
    fprintf('---global:\n')
    for i = 1:4
        [D, B, ~, ~, avgOverlap, DiscErr, ~] = ...
            main(FockDim, Operation, ParamsBounds, Pi, eta, M, C, N, SEQ(i, :), FoM);
        fprintf('%.4f\t%.4f\t%.4f\t%.4f\t[', B, D, DiscErr, avgOverlap);
        fprintf('%.3f ', SEQ(i, :));
        fprintf(']\n');
    end
    fprintf('---local over D:\n');
    [D, B, ~, ~, avgOverlap, DiscErr, ~] = ...
        main(FockDim, Operation, ParamsBounds, Pi, eta, M, C, N, NaN, FoM);
    fprintf('%.4f\t%.4f\t%.4f\t%.4f\n', B, D, DiscErr, avgOverlap);    
    fprintf('\n');
    
end
%}

%% PrintFileInfo
%{
matfile = cellstr(['20141222T151417'; '20141224T204741'; '20141223T185621'; '20141225T024522'; '20141223T201836'; '20141227T140928'; '20141224T180934'; '20141224T204426'; '20150205T140036']);

fprintf('There are %d files.\n\n', length(matfile));

for j = 1:length(matfile)
    PrintFileInfo(['Reporter', matfile{j}]);
end
%}

%% main()
%%{
FockDim = 7;
M = 2; 
C = 2;
N = 2;
eta = 1;
Pi = @APD;
Operation = @CoherentDisplacer;
ParamsBounds = [-1 1];
SEQ = NaN;
FoM_ID = 'D';

[~, ~, ~, ~, ~, ~, ~] = ...
    main(FockDim, Operation, ParamsBounds, Pi, eta, M, C, N, SEQ, FoM_ID);

%}

%% mainGlobal()
%{
FockDim = 3;
M = 2; 
C = 2;
N = 2;
eta = 1;
Pi = @APD;
Operation = @HadamardGate;
ParamsBounds = [-pi pi 3];

[~, ~, ~, ~, mainGlobalClock] = mainGlobal(FockDim, Operation, ...
    ParamsBounds, Pi, eta, M, C, N);
fprintf('mainGlobal was completed in %.2f minutes.', mainGlobalClock/60);
%}

%% Reporter()
%{
C = 2;
FoM_ID = 'D';
%{
M = 2;
FockDim = 3;
minN = 1;
maxN = 10;
etaSamples = 0; % 5. Recall: The actual number of samples is etaSamples+1.
etaMin = 1; % 0.5.
etaMax = 1; % 1.

Operation = @HadamardGate;
Pi = @HD;
ParamsBounds = [-pi/3 pi/3 40];
fprintf('Computing a:\n');
Reporter(FockDim, Operation, ParamsBounds, Pi, M, C, minN, ...
    maxN, etaSamples, etaMin, etaMax, FoM_ID);

etaSamples = 5; 
etaMin = 0.5;
etaMax = 1;
Pi = @APD;
fprintf('Computing b:\n');
Reporter(FockDim, Operation, ParamsBounds, Pi, M, C, minN, ...
    maxN, etaSamples, etaMin, etaMax, FoM_ID);

Pi = @PNRD;
M = 3;
% FAILED
fprintf('Computing c:\n');
Reporter(FockDim, Operation, ParamsBounds, Pi, M, C, minN, ...
    maxN, etaSamples, etaMin, etaMax, FoM_ID);

ParamsBounds = pi/3;
fprintf('Computing d:\n');
Reporter(FockDim, Operation, ParamsBounds, Pi, M, C, minN, ...
    maxN, etaSamples, etaMin, etaMax, FoM_ID);

FockDim = 7;
Operation = @CoherentDisplacer;
etaSamples = 0; 
etaMin = 1;
etaMax = 1;
Pi = @HD;
M = 2;
ParamsBounds = 0.25;
fprintf('Computing e:\n');
Reporter(FockDim, Operation, ParamsBounds, Pi, M, C, minN, ...
    maxN, etaSamples, etaMin, etaMax, FoM_ID);

Pi = @APD;
etaSamples = 5; 
etaMin = 0.5;
etaMax = 1;
fprintf('Computing f:\n');
Reporter(FockDim, Operation, ParamsBounds, Pi, M, C, minN, ...
    maxN, etaSamples, etaMin, etaMax, FoM_ID);

% FAILED at 75%
ParamsBounds = 0;
fprintf('Computing g:\n');
Reporter(FockDim, Operation, ParamsBounds, Pi, M, C, minN, ...
    maxN, etaSamples, etaMin, etaMax, FoM_ID);
%}

FockDim = 7;
Operation = @CoherentDisplacer;
etaSamples = 5; 
etaMin = 0.5;
etaMax = 1;
minN = 1;
maxN = 7;
Pi = @PNRD;
M = 3;
ParamsBounds = 0;
fprintf('Computing h:\n');
Reporter(FockDim, Operation, ParamsBounds, Pi, M, C, minN, ...
    maxN, etaSamples, etaMin, etaMax, FoM_ID);

ParamsBounds = 0.25;
fprintf('Computing i:\n');
Reporter(FockDim, Operation, ParamsBounds, Pi, M, C, minN, ...
    maxN, etaSamples, etaMin, etaMax, FoM_ID);
%}