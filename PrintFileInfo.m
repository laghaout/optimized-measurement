function PrintFileInfo(matfile, N)

matfilePath = ['./reports/', matfile, '.mat'];

if exist(matfilePath, 'file') == 2
    
    load(matfilePath);
    
    if exist('N', 'var')
        minN = min(N);
        maxN = max(N);
    elseif ~exist('minN', 'var')
        minN = 1;
    end
    
    minClock = min(min(mainClock(:, minN:maxN)));
    maxClock = max(max(mainClock(:, minN:maxN)));
    mainClockTotal = sum(sum(mainClock(:, minN:maxN)));
    mainClock = mean(mainClock(:, minN:maxN)); % Mean over all eta's
    b = 2;
    c = mean(log(mainClock)./((minN:maxN)*log(b)));
    fitClock = b.^(c*(minN:maxN))/3600;
    mainClock = mainClock/3600;
    minClock = minClock/3600;
    maxClock = maxClock/3600;

    fprintf('**** %s\n', matfile)
    fprintf('M = %d, C = %d, Operation = %s, Pi = %s\nTotal run time = %.2f hrs, ', ...
        M, C, func2str(Operation), func2str(Pi), mainClockTotal);
    fprintf('eta = [%.1f%%...%.1f%%], N = [%d...%d]\n', ...
        100*etaMin, 100*etaMax, minN, maxN);
    fprintf('ParamsBounds = [');
    fprintf('%f ', ParamsBounds);
    fprintf(']\n');
    fprintf('Average runtime function b^(c*N): b = %f, c = %f\n\n', b, c);
    
else
    
    error('An invalid argument was passed.');
    
end