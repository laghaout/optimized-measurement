%{
Author: Amine Laghaout
Date: 2014-12-09

INPUT

comvars: Cf. main.m
sequence: String of the measurement outcomes
N: Cf. main.m

OUTPUT

sequence: String of the sequence of measurement outcomes written in such a
    way to keep track of its position in the global tree. E.g. in a binary
    tree of depth 3, the sequences '01.' and '100' would represent the
    nodes of indices 5 and 10, respectively, in pre-order traversal. If the
    tree is a subtree, i.e., that's not rooted at the global root, then
    sequences such as '--01.' and '--100' mean that the subtree is of depth
    four, as above, except that it is rooted at level 2. Hence the '--' at
    the beginning indicating the sequence of the first subtree leading up
    to the current subtree.
%}
function sequence = DisplaySequence(comvars, sequence, N)

% The the current tree is a subtree that's rooted at an internal node of 
% the global subtree, then all previously detemined measurement outcomes
% are represented by a dash and all subsequencted undetermined measurement
% outcomes are represented by a dot.
if comvars.RootAt ~= 0
    sequence = [repmat('-', 1, comvars.RootAt), sequence];
    sequence = [sequence, repmat('.', 1, N - length(sequence))];
    
% If the current tree has the same root as the global tree, then all
% subsequent measurement yet to be determined are represented by dots
else
    sequence = [sequence, repmat('.', 1, N - length(sequence))];   

end
