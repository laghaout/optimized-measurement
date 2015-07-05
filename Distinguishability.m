%{
Author: Amine Laghaout
Date: 2014-12-12

INPUT

P: Probabilities of the potential states (rows) at each of the leaves 
    (columns)
p: Array of classical probabilities for each potential state

OUTPUT

D: Distinguishability based on the Hellinger distance (cf. Wikipedia)
%}
function D = Distinguishability(P, p)

C = length(p);  % Number of potential states
WBC = 0;        % Weighted Bhattacharyya coefficient
ImTol = 1e-6;   % Maximum imaginary amplitude tolerated (1e-6 by default)

% Add the Bhattacharyya coefficients pairwise between the PDFs of all the 
% potential states, making sure not to count duplicates
for i = 1:C
    for j = 1:C      
        
        if i ~= j
            WBC = WBC + sqrt(p(i)*P(i,:))*sqrt(p(j)*P(j,:))'/(1 - p(i));
        end
        
    end
end

% Return the Hellinger distance
D = sqrt(1 - WBC);

if D < 0 || D > 1 || isnan(D) || abs(imag(D)) > ImTol
    error('Invalid distinguishability D = %s%% (D = %e + i*%e)', num2str(100*D), real(D), imag(D));
end
