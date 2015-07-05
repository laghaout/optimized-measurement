%{
Author: Amine Laghaout
Date: 2014-12-15

Given a two-mode state 'rhotilde' spanning two-Hilbert spaces of dimensions
'FockDim' each, return the single-mode state 'rho' that's transmitted when
the second mode is traced out.
%}
function rho = PartialTrace(FockDim, rhotilde)

rho = NaN(FockDim);

for n = 1:FockDim
    for m = 1:FockDim
        rho(n, m) = trace(rhotilde((n-1)*FockDim+1:n*FockDim, ...
            (m-1)*FockDim+1:m*FockDim));
    end
end

