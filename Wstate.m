%{
Author: Amine Laghaout
Date: 2014-12-17

This function returns the density matrix of a W state over 'NumModes' of 
dimension 'FockDim' each.

THIS IS ONLY IMPLEMENTED FOR 2, 3, OR 4 MODES.
%}
function rho = Wstate(FockDim, NumModes)

v0 = zeros(FockDim, 1); % Fock vector of the vacuum
v1 = zeros(FockDim, 1); % Fock vector of the single photon
v0(1:length(Fock(0,1))) = Fock(0,1);
v1(1:length(Fock(1,1))) = Fock(1,1);

switch NumModes
    case 2
        v = kron(v1, v0) + kron(v0, v1);
    case 3
        v = kron(kron(v1, v0), v0) + kron(kron(v0, v1), v0) + kron(kron(v0, v0), v1);
    case 4
        v = kron(kron(kron(v1, v0), v0), v0) + kron(kron(kron(v0, v1), v0), v0) + kron(kron(kron(v0, v0), v1), v0) + kron(kron(kron(v0, v0), v0), v1);
    otherwise
        error('Wrong number of modes');
end

rho = v*v';
rho = rho/trace(rho);