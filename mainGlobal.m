%{
Author: Amine Laghaout
Date: 2014-12-27

This function looks for the global optimal set of parameters over all nodes
of the tree where the parameters are sampled from 'ParamsBounds'. All other
arguments and return values are the same as those defined in main.m.

The optimal values found are for the distinguishability 'D', the Bell
factor 'B', the average overlap 'avgOverlap' of the probabilities at the
leaves, and the discrimination error of the superoperators 'DiscErr'. Each
of these may be optimized by different combinations of parameters. The
optimal value and the corresponding combiation of parameters is stored in
a structure, for example, in 'D_opt' and 'D_params', respectively.

One then compares these results with those obtained from the local
optimization strategy of a simple call to main(). (Set 'comvars.VERBOSE' to
true in the latter case to see the set of optimal parameters. An example of
such comparison would be of 
[D, B, avgOverlap, DiscErr, mainGlobalClock] = mainGlobal(7, @CoherentDisplacer, [-1 1 numParams], @APD, 1, 2, 2, 2);
with 
[D, B, ~, ~, avgOverlap, DiscErr, mainClock] = main(7, @CoherentDisplacer, [-1 1 numParams], @APD, 1, 2, 2, 2, NaN);
Remember to change AdaptiveDecision in the call to main() so as to not look
for convergence and so as to not iterate but simply evaluate the array of
parameters.
%}
function [D, B, avgOverlap, DiscErr, mainGlobalClock] = ...
    mainGlobal(FockDim, Operation, ParamsBounds, Pi, eta, M, C, N)

%% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mainGlobalClock = tic;          % Timer
StateDim = NaN(C, 1);           % Fock dimension of the candidate states
TraceErr = 1e-3;                % Maximum tolerated error from unit trace

% Check the validity of the parameter range. The actual number of parameter
% values is 'numParams'+1

switch length(ParamsBounds)

    % The operation Â is to be optimized over a range whose number of
    % parameters is specified by default, namely 10.
    case 2
        minParams = ParamsBounds(1);
        maxParams = ParamsBounds(2);
        numParams = 10;

    % The operation Â is to be optimized over a range whose number of
    % parameters is passed as an argument.
    case 3
        minParams = ParamsBounds(1);
        maxParams = ParamsBounds(2);
        numParams = ParamsBounds(3);
        if numParams < 1
            warning('The minimum number of parameters is 2.');
            numParams = 1;
        end
        
    otherwise
        error('Invalid array of parameter bounds');
end

%% Create the array of parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ParamsArray = (minParams:(maxParams-minParams)/numParams:maxParams)';

% Check that the dimension of the Fock space is large enough to maintain a
% near-unit trace of the density matrices (up to an error 'TraceErr'). This
% is based on the trace of the vacuum once transformed by operation
% 'Operation'.
   
FockDim = SetFockDim(func2str(Operation), ...
    max(abs(minParams), abs(maxParams)), TraceErr, FockDim);

%% Input states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample C equally weighted and equally separated points on a given 
% circumference of the Bloch sphere. In this case we sample along real-
% valued latitudes.

RootState = struct([]);

for i = 1:C % For each potential state...
    
    RootState(i).rho = Qubit(2*pi*(i-1)/C + pi/2, 0); % Fock representation
    RootState(i).P = 1/C;                             % Probabiliy
    StateDim(i) = length(RootState(i).rho);           % Fock dimension
    
end

% Re-adjust the dimension of the Fock space if the states can't fit in.
FockDim = max(max(StateDim), FockDim);

%% Generate the combinations of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lenComb = (M^N-1)/(M-1);            % Length of a parameter combination
numComb = (numParams+1)^lenComb;    % Number of possible such combinations

% Figures of merit and array of parameters to store the final optima
D.opt = -Inf;
B.opt = -Inf;
avgOverlap.opt = Inf;
DiscErr.opt = Inf;
D.params = NaN(1, lenComb);
B.params = NaN(1, lenComb);
avgOverlap.params = NaN(1, lenComb);
DiscErr.params = NaN(1, lenComb);

ParamsMatrix = NaN(numComb, lenComb);

for k = 0:lenComb-1;
    ParamsMatrix(:, k+1) = kron(kron(ones((numParams+1)^k, 1), ...
        ParamsArray), ones((numParams+1)^(lenComb-k-1), 1));
end

%% Distribute the search for the optima
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myPool = gcp('nocreate');
if isempty(myPool)
    myPool = parpool;   
end
S = myPool.NumWorkers-1;
if S > numComb
    S = numComb;
end
numCombPerWorker = floor(numComb/S);

% IMPORTANT: MAKE SURE TO COMMENT OUT THE PARALLELIZATION COMMANDS IN
% main() SO AS TO NOT INTERFERE WITH THE spmd() BELOW.

spmd(S)

    % Figures of merit and array of parameters used by the individual 
    % workers
    
    D_opt = -Inf;
    B_opt = -Inf;
    avgOverlap_opt = Inf;
    DiscErr_opt = Inf;
    D_params = NaN(1, lenComb);
    B_params = NaN(1, lenComb);
    avgOverlap_params = NaN(1, lenComb);
    DiscErr_params = NaN(1, lenComb);

    % Determine the range of rows from 'ParamsMatrix' to be used by the
    % current worker.
    
    mink = (labindex-1)*numCombPerWorker + 1;
    if labindex ~= S
        maxk = mink + numCombPerWorker; % SUBTRACT 1 FROM THIS: TWO WORKERS END UP COMPUTING THIS SAME k
    else
        maxk = numComb;
    end

    for k = mink:maxk

        [currD, currB, ~, ~, curravgOverlap, currDiscErr, ~] = ...
            main(FockDim, Operation, ParamsBounds, Pi, eta, M, C, N, ...
            ParamsMatrix(k, :));
        
        if currD > D_opt
            D_opt = currD;
            D_params = ParamsMatrix(k, :);
        end

        if currB > B_opt
            B_opt = currB;
            B_params = ParamsMatrix(k, :);
        end

        if curravgOverlap < avgOverlap_opt
            avgOverlap_opt = curravgOverlap;
            avgOverlap_params = ParamsMatrix(k, :);
        end

        if currDiscErr < DiscErr_opt
            DiscErr_opt = currDiscErr;
            DiscErr_params = ParamsMatrix(k, :);
        end

        % Display the percentage processed by the current worker
        fprintf('%.1f%% \t', 100*(k-mink)/(maxk-mink));

    end

    fprintf('\n');

end % spmd

% Find the optima over all workers

for k = 1:S
    
    if D_opt{k} > D.opt
        D.opt = D_opt{k};
        D.params = D_params{k};
    end

    if B_opt{k} > B.opt
        B.opt = B_opt{k};
        B.params = B_params{k};
    end    
    
    if avgOverlap_opt{k} < avgOverlap.opt
        avgOverlap.opt = avgOverlap_opt{k};
        avgOverlap.params = avgOverlap_params{k};
    end    

    if DiscErr_opt{k} < DiscErr.opt
        DiscErr.opt = DiscErr_opt{k};
        DiscErr.params = DiscErr_params{k};
    end        
    
end

mainGlobalClock = toc(mainGlobalClock);

ReportID = ['GlobalReporter', datestr(now, 30), '.mat'];
save(ReportID, 'D', 'avgOverlap', 'DiscErr', 'B', 'mainGlobalClock', ...
    'FockDim', 'Operation', 'ParamsBounds', 'Pi', 'eta', 'M', 'C', 'N');