%{
Author: Amine Laghaout
Date: 2014-12-15

INPUT

FockDim: Dimension of the Fock Hilbert space
m: Readout of the PNRD where (m-1) is the number of photons detected
M: Total number of readouts from the PNRD. If m >= M then photon numbers up
    to m-1 can be resolved. Any number of photons equal to m-1 and above 
    cannot be resolved. Setting M = 2 renders this function equivalent to
    APD().
eta: Quantum efficiency. It's a percentage and therefore is represented by
    the power transmission of a fictitious beam splitter.

OUTPUT

Pi = Fock representation of the operation of an photon-number resolving
    detector
%}
function Pi = PNRD(FockDim, m, M, eta)

% Photon numbers up to M - 2 can be determined for sure
if m < M 
    
    Pi = zeros(FockDim);
    k = m-1:FockDim-1;
    Pi(m:FockDim, m:FockDim) = eta^(m-1)*diag((1-eta).^(k-m+1).*arrayfun(@nchoosek,k,(m-1)*ones(1,length(k))));
    
% Photon numbers up equal or greater than M - 1 are indistinguishable.
% I.e., the readout saturates.
elseif m == M

    Pi = eye(FockDim);
    for n = 1:m-1
        Pi = Pi - PNRD(FockDim, n, M, eta);
    end
    
% Invalid readout number
elseif m > M
    
    error('The photon number-resolving detector can only distinguish up to M = %d readouts. Received m = %d.', M, m);
    
end
