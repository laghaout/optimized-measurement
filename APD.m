%{
Author: Amine Laghaout
Date: 2014-12-15

INPUT

FockDim: Dimension of the Fock Hilbert space
m: Readout of the APD (m = 1 for the vacuum, m = 1 for at least one photon)
M: Number of readouts of the APD, which should be M = 2
eta: Quantum efficiency. It's a percentage and therefore is represented by
    the power transmission of a fictitious beam splitter.

OUTPUT

Pi = Fock representation of the operation of an avalanche photo diode
%}
function Pi = APD(FockDim, m, M, eta)

if eta <= 1 && eta >= 0 && M <= 2
    n = (1 - eta).^(0:FockDim-1);
elseif M > 2
    error('Avalanche photo diodes only have two outcomes. Received M = %d.\n', M);
else
    error('Invalid quantum efficiency eta = %s', num2str(eta));
end

switch m
    case 1  % No-click
        Pi = diag(n);
    case 2  % Click
        Pi = eye(FockDim) - diag(n);
    otherwise
       error('Avalanche photo diodes only have two oucomes. Received m = %d.\n', m);
end
