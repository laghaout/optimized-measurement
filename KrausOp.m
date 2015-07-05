%{
Author: Amine Laghaout
Date: 2014-12-17

Apply the 'Kraus' operator on the mode 'Mode' the state 'rho' spanning
'NumModes' states.
%}
function rho = KrausOp(Kraus, rho, Mode, NumModes)

FockDim = length(Kraus);

if Mode < NumModes && Mode > 1
    Kraus = kron(kron(eye(FockDim^(Mode-1)), Kraus), ...
        eye(FockDim^(NumModes-Mode)));
elseif Mode == 1 && NumModes > 1
    Kraus = kron(Kraus, eye(FockDim^(NumModes-1)));
elseif Mode == NumModes && NumModes > 1
    Kraus = kron(eye(FockDim^(NumModes-1)), Kraus);
elseif ~(Mode == NumModes && NumModes == 1)
    error('Invalid number of modes of mode position');
end

rho = Kraus*rho*ctranspose(Kraus);