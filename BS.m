%{
Author: Amine Laghaout
Date: 2014-12-12

INPUT

FockDim: Dimension of the Fock space
t: Amplitude transmission of the beam splitter

OUTPUT

BS: Beam splitter transformation matrix in Fock basis
%}
function BS = BS(FockDim, t)

BS = zeros(FockDim^2);

% WARNING: THIS REQUIRES OPTIMIZATION. USE VECTORIZATION OR arrayfun() AND
% CHECK WHETHER IT ACTUALLY RUNS FASTER. ALSO CONSIDER SAVING THESE BEAM
% SPLITTER MATRICES INSTEAD OF HAVING TO COMPUTE THEM ANEW AT EACH LEVEL.
for n = 0:FockDim-1
    for m = 0:FockDim-1
        for q = 0:n
            for qp = 0:m
                k = q + qp;
                l = n+m-q-qp;
                x = k*FockDim + l + 1;
                y = n*FockDim + m + 1;
                if x <= FockDim^2 && y <= FockDim^2
                    BS(x,y) = BS(x,y) + nchoosek(n,q)*nchoosek(m,qp)*sqrt(factorial(q+qp)*factorial(n+m-q-qp)/(factorial(n)*factorial(m)))*t^(m+q-qp)*(sqrt(1-t^2))^(n-q+qp);
                end
            end
        end
    end
end
