%{
Author: Amine Laghaout
Date: 2014-12-17

Given a certain operation 'Operation' such as coherent displacement or
squeezing of the vacuum, compute the minimum dimension 'FockDim' of the 
Fock space such that the trace of the resulting state's density matrix is
within 'TraceErr' of unity. A preliminary guess of 'FockDim' may be
optionally provided and will be returned unchanged if it is equal or larger
than the minimum Fock space dimension computed. The operation has a certain
parameter 'Param' specific to it.

WARNING: THIS FUNCTION IS INCOMPLETE AND HAS YET TO INCLUDE OPERATIONS
SUCH AS SQUEEZING.
%}
function FockDim = SetFockDim(Operation, Param, TraceErr, FockDim)

% If the minimum Fock space is not specified, start with a dimension of
% two, i.e., large enough to contain a qubit.
if ~exist('FockDim', 'var')
    FockDim = 2;
end

maxn = FockDim-1;

switch Operation
    
    case 'CoherentDisplacer'

        % Compute the trace up to FockDim.
        n = 0:FockDim-1;
        rhoTrace = exp(-abs(Param)^2)* ...
            sum(abs(Param).^(2*n)./arrayfun(@factorial, n));               
        
        % So long as the trace is more than TraceErr away from unity,
        % increment the dimension of the Fock space.
        while 1-rhoTrace > TraceErr
            maxn = maxn + 1;
            rhoTrace = rhoTrace + exp(-abs(Param)^2)* ...
                abs(Param)^(2*maxn)/factorial(maxn);
        end
        
    case 'HadamardGate'
        
        % The Hadamard gate is assumed to only work on qubits. The minimum
        % size of the Fock Hilbert space should therefore be 2. (I.e., the
        % maximum photon number should be 1.)
        maxn = max(1, FockDim - 1);
        
end

FockDim = maxn + 1 + 1;