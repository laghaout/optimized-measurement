%{
Author: Amine Laghaout
Date: 2015-02-08

INPUT

FockDim: Dimension of the Fock Hilbert space     
Operation: Unitary operation used just before the projection
ParamsBounds: Bounds of the parameters in the operation. It consists of an
    array with either 1, 2, or 3 elements. If it consists of a single
    element then the parameter is fixed to that element. If it consists of
    two elements, then they represent the lower and upper bounds of the
    parameters of which tens equally spaced values are sampled. if it
    consists of three parameters, then the first two are the lower and
    upper bound, respectively, and the third specifies the number of values
    sampled.
Pi: Projection operator
eta: Quantum efficiency (percentage) of the measurement device. It is 
    represented by a beam splitter of transmission sqrt(eta).
M: Number of edges per vertex
C: Number of potential states from which the discrimination is based.
N: Number of recursions
SEQ: Combination of the parameters to be applied to the operations in the
    internal nodes. Set to NaN if no such pre-defiend combination of
    parameters is to be used but rather is to computed locally from the
    adaptive decision.
FoM: Figure of merit by which the adaptation is performed.

OUTPUT

D: Distinguishability of the potential states at the leaves.
Pi1: Superoperator corresponding to the first potential state
Pi2: Superoperator corresponding to the second potential state
avgOverlap: Average of the ratio of the minimum to maximum probabilities at
    each leaf
DiscErr: Probabiliy of error resulting of applying Pi1 and Pi2 to the
    input states
mainClock: Timer of the whole function
%}
function [D, B, Pi1, Pi2, avgOverlap, DiscErr, mainClock] = ...
    main(FockDim, Operation, ParamsBounds, Pi, eta, M, C, N, SEQ, FoM_ID)

%% Parallelization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myPool = gcp('nocreate');
if isempty(myPool)
    myPool = parpool; % parpool('DTU_Torque');   
end

% The number of workers is set to one, i.e., to serial computation, if we
% are readily provided with a combination of parameters for the whole tree.
% This latter case happens, for example, when we try to traverse the tree
% for global optimization. The argument 'ParamsBounds' is then modified to
% match this new combination for consistency, although it won't be used 
% other than for expanding the dimension of the Fock space, if necessary.
if isnan(SEQ)
    S = myPool.NumWorkers; 
elseif length(SEQ) == (M^N-1)/(M-1)
    S = 1;
    ParamsBounds = [min(SEQ) max(SEQ) length(SEQ)-1];
else
    error('Invalid sequence of parameters');
end

%% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The level k_d at which the distribution should take place is at N/2. This
% way, the average number of internal nodes to be processed is minimized.
% If there aren't enough workers, use as many of them as possible so as to
% match the number of nodes at the corresponding level.

k_d = ceil(N/2);

if M^k_d > S
    k_d = floor(log(S)/log(M));
end

mainClock = tic;                % Timer
StateDim = NaN(C, 1);           % Fock dimension of the candidate states
NumModes = 1;                   % Total number of modes in the states
Mode = 1;                       % Mode which is measured in the state
TraceErr = 1e-3;                % Maximum tolerated error from unit trace
OutcomeSequence = '';           % String of the sequence of outcomes
comvars.MOVABLE_RANGE = false;  % Move the search range of the parameters? 
comvars.VERBOSE = true;      	% Display text upon execution?
PRINT_SUMMARY = true;          % Display summary?
PLOT_PDF = false;               % Plot the PDFs of the states
comvars.BellModes = 3;          % Number of modes in the Bell state
comvars.RootAt = 0;             % Level at which the current subtree starts
comvars.LeafAt = k_d;           % Level at which the current subtree ends
comvars.SEQ = SEQ;              % Sequence of decision paramters 1:(M^N-1)/(M-1);

% Identifier of the figure of merit to be used:
%   'D': Distinguishability
%   'B': Bell factor
%   'DiscErr': Discrimination error
comvars.FoM_ID = FoM_ID;

% Note: comvars serves to bundle together a bunch of parameters and
% settings to be passed to each node. comvars.RootAt and comvars.LeafAt
% keep track of the levels in the global tree where the current subtree has
% its roots and leaves, respectively. We set them to 0 and k_d,
% respectively, because the first subtree, which is processed in serial, is 
% the first one to be processed. It is rooted at the global root (level 0)
% and goes down to the distribution level k_d.

% If we end up with the distribution level at the root or at the leaves,
% then there is no point in parallelization to begin with.

if k_d == 0 || k_d == N
    comvars.LeafAt = N;
elseif k_d > N || k_d < 0
    error('Invalid distribution level k_d = %d\n', k_d);
end

if comvars.VERBOSE
    fprintf('%d/%d worker(s) distributed at level %d [k_d = %d]\n', ...
        M^(comvars.LeafAt*(comvars.LeafAt ~= N)), S, ...
        comvars.LeafAt*(comvars.LeafAt ~= N), k_d);
end

% Check that the dimension of the Fock space is large enough to maintain a
% near-unit trace of the density matrices (up to an error 'TraceErr'). This
% is based on the trace of the vacuum once transformed by operation
% 'Operation'.

if length(ParamsBounds) == 2 || length(ParamsBounds) == 3
    minParams = min(ParamsBounds(1:2));
    maxParams = max(ParamsBounds(1:2));
elseif length(ParamsBounds) == 1
    minParams = ParamsBounds;
    maxParams = ParamsBounds;
end
    
FockDim = SetFockDim(func2str(Operation), ...
    max(abs(minParams), abs(maxParams)), TraceErr, FockDim);

% Note that the update of 'FockDim' based on operation only makes sense if
% 'comvars.MOVABLE_RANGE' is set to false. Otherwise, the new limits of the
% range, as obtained in AdaptiveDecision() may require larger dimensions of
% the Fock space to still maintain unit trace of the transformed vacuum.

%% Input states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample C equally weighted and equally separated points on a given 
% circumference of the Bloch sphere. In this case we sample along real-
% valued latitudes.

RootState = struct([]);

for i = 1:C % For each potential state...
    
    % (2*i-1)*pi/C
    RootState(i).rho = Qubit((2*i-1)*pi/C, 0);  % Fock representation
    RootState(i).P = 1/C;                       % Probabiliy
    StateDim(i) = length(RootState(i).rho); 	% Fock dimension
    
end

% Re-adjust the dimension of the Fock space if the states can't fit in.
FockDim = max(max(StateDim), FockDim);

% Memorize the root states of the gobal tree. (There will be other root
% states, namely one for each subtree.)
globalRootState = RootState;

%% Recursion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The most basic form the superoperator is that which leaves the states
% unchanged, namely the identity.
comvars.Kraus = eye(FockDim);

% Base tree rooted at the global root (processed in serial)
[D, Leaves] = RecursionTree(0, comvars, FockDim, OutcomeSequence, ...
    globalRootState, Operation, ParamsBounds, Pi, eta, M, C, N);

% Sanity-check the dimensions of Kraus operators at the leaves.
%{
[nrows, ncols, ndepth] = size(Leaves.Kraus);
fprintf('Done with the root tree. Leaves.Kraus is %dx%dx%d.\n', nrows, ncols, ndepth);
%}

% Subtrees (proccessed in parallel)

if k_d > 0 && k_d < N
    
    comvars.RootAt = k_d;   % New roots at the distribution level
    comvars.LeafAt = N;     % New leaves are global leaves
    
    % Save the Kraus operators of the base tree.
    BaseLeavesKraus = Leaves.Kraus;    
    
    % Because of some issues using composites at the end of the spmd(),
    % `comvars' causes an error when attempting to access `comvars.VERBOSE'
    % further down. We shall therefore use a copy `comvars_bis' to be 
    % passed to RecursionTree().
    comvars_bis = comvars; 
    
    spmd(M^k_d) % Distribute the subtrees
        
        % The new root states and their respective probabilities are 
        % recovered from the leaves of the serial subtree.
        for i = 1:C
            RootState(i).rho = Leaves.state(:, :, i, labindex);
            RootState(i).P = Leaves.P(i, labindex);
        end        
        
        % Save the superoperator obtained so far. We'll build upon it.
        comvars_bis.Kraus = BaseLeavesKraus(:, :, labindex);
        
        % Start the recursion for each subtree and save the probabilities
        % and Kraus operators.
        [~, NewLeaves] = RecursionTree(0, comvars_bis, FockDim, ...
            OutcomeSequence, RootState, Operation, ParamsBounds, Pi, ...
            eta, M, C, N);

        % Sanity-check the dimensions of Kraus operators at the leaves.
        %{
        [nrows, ncols, ndepth] = size(NewLeaves.Kraus);
        fprintf('Done with the sub-tree number %d. Leaves.Kraus is %dx%dx%d.\n', labindex, nrows, ncols, ndepth);        
        %}
        
        % Save the new leaves.
        NewLeavesP = NewLeaves.P;
        NewLeavesKraus = NewLeaves.Kraus;
        
    end % spmd()
    
    % Reuse Leaves for the global tree
    clear Leaves;       
    Leaves.P = [];
    Leaves.Kraus = [];    
    
    % Stitch together the Kraus operators and the probabilities at the
    % leaves (i.e., "horizontally" from the different workers).
    for leaf = 1:M^k_d

        Leaves.P = cat(2, Leaves.P, NewLeavesP{leaf});
        Leaves.Kraus = cat(3, Leaves.Kraus, NewLeavesKraus{leaf});

    end
    
    D = Distinguishability(Leaves.P, [globalRootState.P]);
    
end % Subtrees

%% Superoperators for discrimination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WARNING: THIS SECTION IS HARD-CODED FOR TWO SUPEROPERATORS ONLY, LABELLED
% Pi1 and Pi2. IT SHOULD BE EVENTUALLY ADAPTED TO PRODUCE AS MANY
% SUPEROPERATORS AS THERE ARE STATES, I.E. C.

leafPoints = 1:M^N;         % Indices of the leaves of the overall tree
LeafOutcome = NaN(1, M^N);  % Indices of states most likely at each leaf

% Average ratio of min-to-max probability at each leaf. It should ideally
% be 0; in the worst case it is 1, when min = max.
avgOverlap = 0;             

Pi1 = [];   % Superoperator corresponding to RootState(1).rho
Pi2 = [];   % Superoperator corresponding to RootState(2).rho

% A given path down the decision tree is represented by sum of nested Kraus
% operators acting on the initial state. This gives a superoperator for
% that particular path. The overall superoperators Pi1 and Pi2 are
% essentially assembled from these path superoperators for each leaf
% corresponding to either RootState(1) or RootState(2), respectively. Refer
% to SuperOp().

for leaf = leafPoints % At each leaf...
        
    % Idenitify which candidate state produces the maximum probability and
    % which produces the minimum probability.
    [minProb, ~] = min(Leaves.P(1:C, leaf));
    [maxProb, MostLikelyState] = max(Leaves.P(1:C, leaf));
    LeafOutcome(leaf) = MostLikelyState;
    
    % Incorporate the appropriate set of Kraus operators to either the Pi1
    % or Pi2 superoperator if the probabilities are different. If the
    % probabilities are equal, assign them to Pi1 by default but print a
    % warning.
    if abs(maxProb-minProb) > 1e-15
        switch MostLikelyState
            case 1
                Pi1 = cat(3, Pi1, Leaves.Kraus(:,:,leaf));
            case 2
                Pi2 = cat(3, Pi2, Leaves.Kraus(:,:,leaf));
            otherwise
                error('Invalid index of the most likely state at leaf %d', leaf);
        end
    else
        if comvars.VERBOSE || PRINT_SUMMARY
            warning('Equal probabilities at leaf %d (N = %d)\n', leaf, N);
        end
        
        Pi1 = cat(3, Pi1, Leaves.Kraus(:,:,leaf));        
        
        % Even if we assign the superoperator to `Pi1', indicate that there
        % was a tie in the probabilities with a zero (i.e., inconclusive)
        % outcome.
        LeafOutcome(leaf) = 0;
    end 
    
    % The average overlap is computed as the weighted sum of ratio of
    % minimum to maximum.
    avgOverlap = avgOverlap + sum(Leaves.P(:, leaf))*minProb/maxProb;

end % Loop over leaves

% Compute the error probability

rho1 = zeros(FockDim);
rho2 = zeros(FockDim);

rho1(1:StateDim(1), 1:StateDim(1)) = globalRootState(1).rho;
rho2(1:StateDim(2), 1:StateDim(2)) = globalRootState(2).rho;

Pi1rho1 = trace(SuperOp(Pi1, rho1, Mode, NumModes));
Pi1rho2 = trace(SuperOp(Pi1, rho2, Mode, NumModes));
Pi2rho1 = trace(SuperOp(Pi2, rho1, Mode, NumModes));
Pi2rho2 = trace(SuperOp(Pi2, rho2, Mode, NumModes));

DiscErr = globalRootState(1).P*Pi2rho1 + globalRootState(2).P*Pi1rho2;

% Bell factor with the tri-partite W state
B = BellTest(FockDim, comvars.BellModes, func2str(Pi), M, eta, Pi1, Pi2);

%% Determine the optimal assignment of the leaves to the superoperator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
% NOTE THAT THIS CANNOT BE RUN WITH batch().

delete(gcp('nocreate'));
[optB, LeafOutcomeB, optDiscErr, LeafOutcomeDiscErr] = ...
    Leaf2SuperOp(FockDim, comvars, Pi, M, N, eta, Leaves.Kraus, ...
    globalRootState, StateDim, Mode, NumModes);

% THERE WILL BE AN ERROR IN WHAT IS BELOW BECAUSE OF SOME COMPOSITES

if comvars.VERBOSE
    fprintf('**** Optimization of leaf assignments\n');

    if B == optB && DiscErr == optDiscErr
        fprintf('** All optima concur!\n');
    else
        if B ~= optB
            fprintf('** Mismatch in Bell factor\n');
        end
        if DiscErr ~= optDiscErr
            fprintf('** Mismatch in discrimination error\n');
        end
    end
    
    fprintf('Optimal B: %f\nLeaf assignment:\t', optB)
    fprintf('%d\t', LeafOutcomeB(:));
    fprintf('\n');
    fprintf('Optimal DiscErr: %f\nLeaf assignment:\t', optDiscErr)
    fprintf('%d\t', LeafOutcomeDiscErr(:));
    fprintf('\n');    
    
    fprintf('****\n');    
end
%}

%% Reports
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mainClock = toc(mainClock);

% Check that the probabilities at the leaves add up to unity
if abs(1-sum(sum(Leaves.P))) > 1e-2
    error('Total probability %f%% at the leaves is not unity.\n', ...
        'Consider increasing the dimension of the Hilbert space.', ...
        100*sum(sum(Leaves.P)));

% Check that the superoperator probabilities add up to unity
elseif abs(Pi1rho1+Pi2rho1 - 1)  + abs(Pi1rho2+Pi2rho2 - 1) > 1e-2
    error('Superoperator probabilities do not add up to unity: Pi1rho1+Pi2rho1 = %.2f%% and Pi1rho2+Pi2rho2 = %.2f%%', 100*(Pi1rho1+Pi2rho1), 100*(Pi1rho2+Pi2rho2));    
    
% If all probabilities figure, print the summary
elseif PRINT_SUMMARY
    if comvars.VERBOSE
        fprintf('Leaves:\t\t');
        fprintf('%d\t', 1:M^N);
        fprintf('\nOutcome:\t');
        fprintf('%d\t', LeafOutcome(:));
        fprintf('\n');
        fprintf('Pi1rho1: %.2f\tPi1rho2: %.2f\nPi2rho1: %.2f\tPi2rho2: %.2f\n', 100*Pi1rho1, 100*Pi1rho2, 100*Pi2rho1, 100*Pi2rho2);        
    end
    fprintf('FockDim = %d, M = %d edges per vertex, C = %d states, N = %d recursions, FoM = %s\n', FockDim, M, C, N, comvars.FoM_ID);
    fprintf('D = %.3f%% B = %.5f avgOverlap = %.3f%% DiscErr = %.3f%%\n', 100*D, B, 100*avgOverlap, 100*DiscErr);
    fprintf('Run time = %.2f sec\n', mainClock);
    
end

%% Probabiliy distribution functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PLOT_PDF
    
    close all;
    set(0, 'DefaultFigureWindowStyle', 'docked');

    bar(leafPoints, 100*Leaves.P', 'BarLayout', 'stacked', 'BarWidth', 1);
    colormap copper;
    xlabel('Leaf');
    ylabel('Probability [%]');
    axis tight;
    title({['D = ', num2str(100*Distinguishability(Leaves.P, ...
        [globalRootState.P]), '%2.1f'), '%'], ['M = ', num2str(M), ...
        ', C = ', num2str(C), ', N = ', num2str(N)]});
    print('-dpng', ['.\figures\PDF_', func2str(Operation), '_', ...
        func2str(Pi),'_M_', num2str(M), '_C_', num2str(C), '_N_', ...
        num2str(N), '.png']);
    L = cell(1, C);
    for i = 1:C
        L(1, i) = cellstr(['\rho_', num2str(i)]);
    end
    legend(L);

end