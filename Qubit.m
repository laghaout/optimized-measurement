%{
Author: Amine Laghaout
Date: 2014-12-15

Return the normalized qubit cos(theta/2)|0> + sin(theta/2)*exp(1i*phi)|1>
%}
function rho = Qubit(theta, phi)

rho = [cos(theta/2)^2 exp(-1i*phi)*cos(theta/2)*sin(theta/2); ...
    exp(1i*phi)*sin(theta/2)*cos(theta/2) sin(theta/2)^2];
