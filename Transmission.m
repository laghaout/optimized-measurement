%{
Author: Amine Laghaout
Date: 2014-12-12

INPUT

k: Level in the global tree
N: Depth of the global tree

OUTPUT

t: Amplitude ransmission of the beam splitter such that the input state is 
    equally split over all N modes
%}
function t = Transmission(k, N)

if k <= N
    
    t = sqrt((N - k)/(N - k + 1));
    
    if abs(t) > 1 || isnan(t) || imag(t) ~= 0
        error('Invalid transmission t = %s at k = %d for N = %d', num2str(t), k, N);
    end
    
else
    
    error('Level k = %d > depth N = %d', k, N);

end
