%{
Author: Amine Laghaout
Date: 2014-12-15

Hadamard operator on a qubit expressed in a Hilbert space of dimension 
'FockDim'. The angle of rotation is 'theta'.
%}
function H = HadamardGate(FockDim, theta)

H = eye(FockDim);
H(1:2, 1:2) = [cos(theta/2) -sin(theta/2); sin(theta/2) cos(theta/2)];