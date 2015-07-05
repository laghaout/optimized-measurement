%{
Author: Amine Laghaout
Date: 2014-12-15

Consider two Fock numbers [n m] stored in `braket' and an operation
`UnitaryOperation'. This function returns the reduced matrix 
<n|UnitaryOperation|m> as `ReducedOperation' where the bra and ket pertain
to the second mode.

To render this functiom more general, the reduction <n|UnitaryOperation|m>
may actually act over a set of n's and m's, each with a given weight stored
in an array `weight' with the same number of rows as 'braket'. For example:
braket = [1 0; 5 2] and weight = [1 2]' means that we return 
(1/3)*<1|UnitaryOperation|0> + (2/3)*<5|UnitaryOperation|2>
%}
function ReducedOperation = ProjectMode2(braket, weight, UnitaryOperation)

% For our purposes, `UnitaryOperation' is assumed to span two modes of
% Hilbert space dimension FockDim each. 
FockDim = sqrt(length(UnitaryOperation));

% Create a matrix A to act as a stencil for the bras and kets projected
% upon in the second mode.
A = zeros(FockDim);
weightSum = sum(weight);
[braketRows, ~] = size(braket);
for i = 1:braketRows
    A(braket(i,1)+1, braket(i,2)+1) = weight(i, 1)/weightSum;
end

% Expand to the full two-mode extent and "punch out" all elements in
% UnitaryOperation that don't match the bras and kets of braket.
UnitaryOperation = UnitaryOperation.*repmat(A, FockDim, FockDim);

% Trace out; i.e., collapse the submatrices of the second mode.
ReducedOperation = NaN(FockDim);
for n = 1:FockDim
    for m = 1:FockDim
        ReducedOperation(n, m) = ...
            sum(sum(UnitaryOperation((n-1)*FockDim+1:n*FockDim, ...
            (m-1)*FockDim+1:m*FockDim)));
    end
end

