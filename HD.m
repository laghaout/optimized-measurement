%{
Author: Amine Laghaout
Date: 2014-12-15

INPUT

FockDim: Dimension of the Fock Hilbert space
m: Readout of the homodyne detector (HD). In the present case m = 1 and m =
    2 represent negative and positive quadratures respectively.
M: Total number of readouts from the HD. 
eta: Quantum efficiency. It's a percentage and therefore is represented by
    the power transmission of a fictitious beam splitter. THE CURRENT
    VERSION CAN ONLY ACCOMMODATE UNIT QUANTUM EFFICIENCIES.

OUTPUT

Pi = Fock representation of the operation of homodyne detector that
    bins between positive and negative quadratures.
%}
function Pi = HD(FockDim, m, M, eta)

if eta ~= 1
    error('The HD has not been implemented for non-unit quantum efficiencies yet.')
elseif M > 2
    error('The HD is implemented to only have two quadrature intervals: x < 0 and x >0. Received M = %d.\n', M);
end

switch m
    case 1
        Pi = csvread('./HDnegative.csv', 0, 0, [0 0 FockDim-1 FockDim-1]);
    case 2
        Pi = csvread('./HDpositive.csv', 0, 0, [0 0 FockDim-1 FockDim-1]);
    otherwise
        error('The outcome index m = %d for a binary homodyne detector is invalid.', m);
end
