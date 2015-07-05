%{
Author: Amine Laghaout
Date: 2015-02-08

This function, given a complete set of Kraus operators from the decision
tree, will shuffles the assignments of the leaves to the superoperators of
the candidate qubits in all possible combinations in order to maximim the
Bell factor and the discrimination error.

WARNING: THIS IS HARD-CODED FOR TWO CANDIDATE STATES ONLY.
%}
function [opt_B, LeafOutcomeB, opt_DiscErr, LeafOutcomeDiscErr] = ...
    Leaf2SuperOp(FockDim, comvars, Pi, M, N, eta, LeavesKraus, ...
    globalRootState, StateDim, Mode, NumModes)

leafPoints = 1:M^N;             % Indices of the leaves of the overall tree
numComb = 2^(M^N);              % Number of combinations
Combs = dec2bin(0:numComb-1);   % List of combinations
opt_B = -Inf;                    % Optimal Bell factor
opt_DiscErr = Inf;               % Optimal discrimination error

% Initial qubits at the root
rho1 = zeros(FockDim);          
rho2 = zeros(FockDim);
rho1(1:StateDim(1), 1:StateDim(1)) = globalRootState(1).rho;
rho2(1:StateDim(2), 1:StateDim(2)) = globalRootState(2).rho;

% Array of assignments at each leaf
LeafOutcome = NaN(1, M^N);          % Current assignment
LeafOutcome_B = NaN(1, M^N);        % Assignment optimizing B
LeafOutcome_DiscErr = NaN(1, M^N);	% Assignment optimizing DiscErr

%% Determine the number of combinations assigned to each worker

myPool = gcp('nocreate');
if isempty(myPool)
    myPool = parpool('DTU_Torque', 32);   
end

S = myPool.NumWorkers-1; % Take into account the reference worker

% If there is more workers than combinations, assign a combination to
% each worker
if S > numComb
    S = numComb;
end

numCombPerWorker = floor(numComb/S);

fprintf('Number of workers: %d\nNumber of combinations per worker: %d\n', S, numCombPerWorker);

spmd(S)
    
    % Determine the index of the first combination for this worker
    mink = (labindex-1)*numCombPerWorker + 1;
    if labindex ~= S
        maxk = mink + numCombPerWorker - 1;
    else
        maxk = numComb;
    end
    
    for k = mink:maxk % Range of combinations assigned to this worker
        
        Pi1 = [];   % Superoperator corresponding to RootState(1).rho
        Pi2 = [];   % Superoperator corresponding to RootState(2).rho
        
        for leaf = leafPoints
            
            switch Combs(k, leaf)
                case '0' % Assign to Pi1
                    LeafOutcome(leaf) = 1;
                    Pi1 = cat(3, Pi1, LeavesKraus(:,:,leaf));
                case '1' % Assign to Pi2
                    LeafOutcome(leaf) = 2;
                    Pi2 = cat(3, Pi2, LeavesKraus(:,:,leaf));
                otherwise
                    error('Wrong assignment of the superoperator at leaf %d', leaf);
            end
            
            
        end % Loop over the leaves

        % Determine the optimal combinations leading to the highest Bell factor
        % and lowest discrimination error

        Pi1rho2 = trace(SuperOp(Pi1, rho2, Mode, NumModes));
        Pi2rho1 = trace(SuperOp(Pi2, rho1, Mode, NumModes));
        DiscErr = globalRootState(1).P*Pi2rho1 + globalRootState(2).P*Pi1rho2;

        if DiscErr < opt_DiscErr
            opt_DiscErr = DiscErr;
            LeafOutcome_DiscErr = LeafOutcome;
        end    

        B = BellTest(FockDim, comvars.BellModes, func2str(Pi), M, eta, Pi1, Pi2);    

        if B > opt_B
            opt_B = B;
            LeafOutcome_B = LeafOutcome;
        end
        
        %{
        fprintf('%f\t%f\t%f\t%f:\t', B, opt_B, DiscErr, opt_DiscErr);
        fprintf('%d \t', LeafOutcome);
        fprintf('\n');
        %}
        
    end % Loop over combinations
        
end % spmd()

% Determine the optima in Bell factor and discrimination error over all
% workers

clear optB optDiscErr LeafOutcomeB LeafOutcomeDiscErr
optB = -Inf;
optDiscErr = Inf;
LeafOutcomeB = NaN(1, M^N); 
LeafOutcomeDiscErr = NaN(1, M^N); 

for w = 1:S % For all workers...
    
    if opt_B{w} > optB
        optB = opt_B{w};
        LeafOutcomeB = LeafOutcome_B;
    end

    if opt_DiscErr{w} < optDiscErr
        optDiscErr = opt_DiscErr{w};
        LeafOutcomeDiscErr = LeafOutcome_DiscErr;
    end
        
end

fprintf('optB = %.5f\noptDiscErr = %.2f%%\n', optB, 100*optDiscErr);
