%{
Author: Amine Laghaout
Date: 2015-01-30

INPUT

NodeID: Index of the root node in pre-order traversal
comvars: Cf. main.m
FockDim: Cf. main.m
OutcomeSequence: Cf. main.m
RootState: Density matrix and classical probabilities of the potential
    states input at the root node
Operation: Cf. main.m
ParamBounds: Cf. main.m
Pi: Cf. main.m
eta: Cf. main.m
M: Cf. main.m
C: Cf. main.m
N: Cf. main.m

OUTPUT

D: Distinguishability of the potential states input at the root via
    RootState once they reach the leaves
Leaves: Structure containing the information about the leaves, namely,
    their index, the probability of the potential states reaching them, the
    density matrix of those states, and the set of Kraus operators
    transforming the states from the root to any given leaf
%}
function [D, Leaves] = RecursionTree(NodeID, comvars, FockDim, ...
    OutcomeSequence, RootState, Operation, ParamsBounds, Pi, eta, M, C, N)

% Depth of the current tree
currN = comvars.LeafAt - comvars.RootAt;    

% Probabilities of the C potential states at the M^currN leaves
Leaves.P = NaN(C, M^currN);

% Index of the leaves. Starts with 1 and ends with M^currN
Leaves.ID = 1;

% Fock representation of the states as FockDim-by-FockDim matrices. There
% are C such matrices for each potential state at ech M^currN leaf.
Leaves.state = zeros(FockDim, FockDim, C, M^currN);

% Set of Kraus operators leading up to the leaves. There are currN such
% nested matrices of dimension FockDim. The horizontal dimension FockDim^2
% is the result of tracing out the mode at each level. There are as many
% sets of Kraus operations as there are leaves, namely M^currN. Note
% however, that if if we are dealing with subtree, they need to carry the
% whole sequence of Kraus operators including those belonging the the root
% tree. In that case, we need to allocate `N' levels of Kraus operators,
% not just `currN'.
if comvars.RootAt == 0
    Leaves.Kraus = NaN(currN*FockDim, FockDim^2, M^currN);
else
    Leaves.Kraus = NaN(N*FockDim, FockDim^2, M^currN);
end

if comvars.VERBOSE
    fprintf('ID\tk\tsequence\tP\t');
    fprintf('P%d\t', 1:C);
    fprintf('FoM\t\tparam\t\tremarks\n');
end

% Start the recursion from the root
[~, Leaves] = RecursionNode(comvars, Leaves, FockDim, NodeID, ...
    OutcomeSequence, N, M, comvars.RootAt, RootState, eta, Operation, ...
    ParamsBounds, Pi);

% Compute the distinguishability of the potential states at the leaves
% based on initial classical probabilities RootState.P
D = Distinguishability(Leaves.P, [RootState.P]);

% If we started at the global root, then we should make sure that the
% probabilities add up to unity.
if comvars.RootAt == 0 && abs(1-sum(sum(Leaves.P))) > 1e-2
    error('Total probability %f%% at the leaves is not unity.\nConsider increasing the dimension of the Hilbert space.', 100*sum(sum(Leaves.P)));
elseif comvars.VERBOSE
    fprintf('Total probability at the leaves = %f%%\n', 100*sum(sum(Leaves.P)));
end
