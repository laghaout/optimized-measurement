%{
Author: Amine Laghaout
Date: 2014-12-17

This function returns the Fock vector of a pure superposition of photon
number states stored in 'photonNums' with respective coefficients 'coeffs'.
The coefficients are scaled so as to return a normalized density matrix.
%}
function v = Fock(photonNums, coeffs)

if length(photonNums) == length(coeffs)
    
    v = zeros(length(coeffs), 1);
    v(photonNums+1) = coeffs/sqrt(sum(abs(coeffs).^2));
    
else
    
    error('Mismatch between the photon numbers and their coefficients');
    
end