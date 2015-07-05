%{
Author: Amine Laghaout
Date: 2014-12-17

Compute the Bell factor for a W state as defined in Phys. Rev. A 84, 062127
(2011) [doi:10.1103/physreva.84.062127].

There are two types of projections: Z and X. The Z projection is always
performed with a photon counting device such as an APD. The X projection is
implemented either via the method 'PiMethod' or using superoperators 'Xm'
and 'Xp', corresponding to the negative and positive projections onto X,
respectively. The number of modes in the W state is 'NumModes', whereby
each mode is represented by in a Fock space of dimension 'FockDim'. 'M', as
defined in main.m, is the number of possible readouts of the measurement
device corresponding to 'PiMethod'.
%}
function BellFactor = BellTest(FockDim, NumModes, PiMethod, M, eta, Xm, Xp)

%% If neither superoperator Xm or Xp is provided, use pre-defined ones

if ~exist('Xp', 'var') || ~exist('Xm', 'var')

    switch upper(PiMethod)
        case 'IDEAL'
            Xp = zeros(FockDim);
            Xm = zeros(FockDim);
            Xp(1:2, 1:2) = Qubit(pi/2,0);
            Xm(1:2, 1:2) = Qubit(-pi/2,0);
        case 'HD'
            Xp = HD(FockDim, 2, M, eta); 
            Xm = HD(FockDim, 1, M, eta);
        case 'APD'
            Xp = APD(FockDim, 2, M, eta); 
            Xm = APD(FockDim, 1, M, eta);
        case 'CD_APD'
            error('Not implemented yet');
            % Load the stored superoperators for the given N, M, eta,
            % PiMethod...
        case 'HG_APD'
            error('Not implemented yet');
        case 'CD_HD'    
            error('Not implemented yet');
        case 'HG_HD'
            error('Not implemented yet'); 
        otherwise
            error('Invalid measurement method');
    end   

end

if NumModes ~= 2 && NumModes ~= 3 && NumModes ~= 4
    error('Bell test not implemented for %d modes', NumModes);
end

%% Add up the probabilities making up the Bell factor

Zp = APD(FockDim, 1, 2, eta); 
Zm = APD(FockDim, 2, 2, eta);

rhoin = Wstate(FockDim, NumModes);

rho = rhoin;
rho = SuperOp(Zm, rho, 1, NumModes);
for Mode = 2:NumModes
    rho = SuperOp(Zp, rho, Mode, NumModes);
end
P1 = trace(rho);

rho = rhoin;
rho = SuperOp(Xp, rho, 1, NumModes);
rho = SuperOp(Xm, rho, 2, NumModes);
for Mode = 3:NumModes
    rho = SuperOp(Zp, rho, Mode, NumModes);
end
P2 = trace(rho);

rho = rhoin;
for Mode = 1:NumModes
    rho = SuperOp(Xp, rho, Mode, NumModes);
end
P3 = trace(rho);

rho = rhoin;
for Mode = 1:NumModes
    rho = SuperOp(Xm, rho, Mode, NumModes);
end
P4 = trace(rho);

BellFactor = NumModes*P1 - NumModes*(NumModes-1)*P2 - P3 - P4;