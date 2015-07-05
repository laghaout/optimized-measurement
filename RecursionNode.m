%{
Author: Amine Laghaout
Date: 2105-01-30

INPUT

comvars: Cf. main.m
Leaves: Cf. RecursionTree.m
FockDim: Cf. main.m
NodeID: Cf. RecursionTree.m
OutcomeSequence: Cf. main.m
N: Cf. main.m
M: Cf. main.m
k: Level of the node in the overall tree
state: Input state
eta: Cf. main.m
Operation: Cf. main.m
ParamsBounds: Cf. main.m
Pi: Cf. main.m

OUTPUT
 
NodeID: Highest node index of among the children of the current node in
    pre-order traversal
Leaves: Cf. RecursionTree.m
%}
function [NodeID, Leaves] = RecursionNode(comvars, Leaves, FockDim, ...
    NodeID, OutcomeSequence, N, M, k, state, eta, Operation, ...
    ParamsBounds, Pi)

% Recover  the number of potential states
C = length(state);

% Print the pre-order traversal index of the node, the total probability of
% reaching that node given the classical mixture of potential states and
% the probability of each state reaching that node.
if comvars.VERBOSE
    fprintf('%d\t%d\t[%s]\t\t', NodeID, k, ...
        DisplaySequence(comvars, OutcomeSequence, N));
    fprintf('%.2f\t', 100*sum([state.P]));
    fprintf('%.2f\t', 100*[state(1:C).P]);
end

%% Sanity-check that the superoperator reproduces the propagated states.
%{
% Initial state. (Let's look at the first candidate qubit in this case.)
statei = 1; % Use 2 instead if we want to work the the second candidate.
initQubit = zeros(FockDim);
initQubit(1:2, 1:2) = Qubit(2*pi*(statei-1)/C + pi/2, 0); 
initQubitP = 1/2;
UFlow = 1e-15;

if isfield(comvars, 'Kraus')
    
    % Current qubit (i.e. the one that was propagated)
    currQubit = zeros(FockDim);
    currQubit(1:length(state(statei).rho), 1:length(state(statei).rho)) = state(statei).rho;
    currQubitP = state(statei).P;
    currQubit = currQubit/currQubitP;
    
    % Qubit produced by the superoperation
    KrausQubit = SuperOp(comvars.Kraus, initQubit, 1, 1);
    KrausQubitP = trace(KrausQubit)*initQubitP;
    KrausQubit = KrausQubit/KrausQubitP;

    if sum(sum(KrausQubit-currQubit)) > UFlow
        error('Evolved qubit does not match the one that is Kraus-operated: %f', sum(sum(KrausQubit-currQubit)));
    elseif abs(currQubitP-KrausQubitP) > UFlow
        error('Mismatch of probabilities: %.2f%% vs %.2f%%', 100*currQubitP, 100*KrausQubitP);
    else
        [KopRows, KopCols] = size(comvars.Kraus);
        fprintf('[Kop=%dx%d currP=%.2f%% KrausP=%.2f%%] ', ...
            KopRows, KopCols, 100*currQubitP, 100*KrausQubitP);            
    end

else

    error('Leaves.Kraus is undefined.');

end
%}

%% Internal node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if k < comvars.LeafAt
    
    % Fock representation of the beam splitting operation
    B = BS(FockDim, Transmission(k+1, N));
        
    % Find the best operation OptimalA to be applied in the measured mode
    % OptimalParam: Optimal parameter of the operation
    % optFoM: Optimized figure of merit (e.g. distinguishability, Bell 
    %   factor)
    % Remarks: Any remarks regarding the performance of the decision
    %   process
    % OptimalA: Fock representation of the optimal operation

    % If we are already provided with a sequence of parameters, just pass
    % the parameter that corresponds to the index of the internal node to
    % the adaptive decision, which will simply propagate the state with
    % that parameter.
    if isnan(comvars.SEQ)
        OptimalParam = ParamsBounds;
    else
        % The index of this internal node, among internal nodes, is
        % NodeID+2-Leaves.ID
        OptimalParam = comvars.SEQ(NodeID+2-Leaves.ID); 
    end

    [OptimalParam, optFoM, Remarks, OptimalA] = AdaptiveDecision(comvars, ...
        FockDim, M, state, B, eta, Operation, Pi, OptimalParam, k);
    
    if comvars.VERBOSE
        fprintf('%f\t%f\t%s\n', optFoM, OptimalParam, Remarks);
    end

    % Proceed with pre-order recursion through each child node emerging
    % from a measurement outcome `m'
    for m = 1:M
        
        % Projection operation corresponding to the outcome `m'
        mPi = Pi(FockDim, m, M, eta);
        m_state =  UpdateNode(state, B, OptimalA, mPi);
        [sqrtPi, ~] = sqrtm(mPi);

        % Overall two-mode operation on the state
        BlackBox = kron(eye(FockDim), sqrtPi)*kron(eye(FockDim), OptimalA)*B;
        
        % Trace out the second mode to obtain the Kraus operator
        for n = 0:FockDim-1
            comvars.Kraus(k*FockDim+1:(k+1)*FockDim, n*FockDim+1:(n+1)*FockDim) = ...
                ProjectMode2([n 0], 1, BlackBox);
        end
        
        % Recursive step
        [NodeID, Leaves] = RecursionNode(comvars, Leaves, FockDim, ...
            NodeID+1, [OutcomeSequence, num2str(m-1)], N, M, k+1, ...
            m_state, eta, Operation, ParamsBounds, Pi);
        
    end 

%% Leaf node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif k == comvars.LeafAt
    
    if comvars.VERBOSE                
        fprintf(' -< %d\n', Leaves.ID);
    end
    
    % Probabilities of the potential states at the leaf
    Leaves.P(:, Leaves.ID) = [state(1:C).P]';   

    % Completed set of Kraus operators
    Leaves.Kraus(:, :, Leaves.ID) = comvars.Kraus;
        
    % Final representation of the potential states
    for i = 1:C
        rhoLength = length(state(i).rho);
        Leaves.state(1:rhoLength, 1:rhoLength, i, Leaves.ID) = state(i).rho;
    end
    
    Leaves.ID = Leaves.ID + 1;

% Error in the level of the node    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

else
    
    error('Negative level k = %d of the tree\n', k);
    
end
