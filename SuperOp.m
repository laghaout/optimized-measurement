%{
Author: Amine Laghaout
Date: 2014-12-17

Function for the superoperation, i.e., the application of the superoperator
'Kraus' on mode 'Mode' of a state 'rhoin' spanning 'NumModes' modes. It
returns the state 'rhoout' that is reduced from 'rhoin'.

Note: By Kraus operator, we really mean the set of Kraus operators making
up the superoperator.

WARNING: Be extra careful about the fact that `FockDim' is reconstituted,
not passed as an argument. This causes issues in some cases.
%}
function rhoout = SuperOp(Kraus, rhoin, Mode, NumModes)

% Recover the dimension of the Hilbert space in any given mode
FockDim = round(length(rhoin)^(1/NumModes));

% Recover the number of recursions from the number of elements in the
% three-dimensional Kraus operator
[row, col, dep]  = size(Kraus);
N = row/(col/FockDim);

%% This is a simple projection: 'Kraus' is a simple projection matrix.

if row == FockDim && col == row

    % 'Kraus' is really just a projection matrix. To get the Kraus operator
    % proper, take its square root.
    [muKraus, ~] = sqrtm(Kraus);
    
    % Apply the Kraus operator
    rhoout = KrausOp(muKraus, rhoin, Mode, NumModes);  
    
%% This is a superoperation: Nested operation of the Kraus operators

elseif (row + col) ~= 0

    rhoin_init = rhoin;
    rhoout = zeros(FockDim^NumModes);
    
    for l = 1:dep % For all leaves...

        rhoin = rhoin_init;

        for k = 1:N % Down all levels...

            rhoTemp = zeros(FockDim^NumModes);               

            for n = 1:FockDim % All basis elements (i.e., trace)...

                % Actual Kraus operator corresponding to the current leaf,
                % level, and basis element
                muKraus = Kraus((k-1)*FockDim+1:k*FockDim, ...
                    (n-1)*FockDim+1:n*FockDim, l);

                % Add up the application of the Kraus operators
                rhoTemp = rhoTemp + ...
                    KrausOp(muKraus, rhoin, Mode, NumModes);

            end

            % Move on to the next level (recursion)
            rhoin = rhoTemp;

        end

        % The output is the sum of outputs from every leaf
        rhoout = rhoout + rhoin;

    end

%% The Kraus operator is empty 
    
else
    
    rhoout = 0;
    
end

