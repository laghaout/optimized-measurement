%{
Author: Amine Laghaout
Date: 2014-12-15

Displacement operator in Fock basis of dimension 'FockDim' and displacement
amplitude 'alpha'.
%}
function D = CoherentDisplacer(FockDim, alpha)

D = NaN(FockDim);

for n = 0:FockDim-1
    for m = 0:FockDim-1
        k = (0:min([m n]))';
        D(n+1,m+1) = sum(sqrt(factorial(n).*factorial(m)).*(-conj(alpha)).^(m-k).*alpha.^(n-k)./(factorial(k).*factorial(m-k).*factorial(n-k)));
    end
end

D = exp(-abs(alpha)^2/2)*D;