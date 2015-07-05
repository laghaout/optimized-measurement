function Tester(alpha, SEQ, x, y, cumP)

N = length(SEQ);

for k = 1:N
    [x, y, P] = Analytical(alpha(k), (Transmission(k, N))^2, x, y, SEQ(k));
    cumP = cumP*P;
    fprintf('Cumulative P = %f%%\n', 100*cumP);
end