%{
Author: Amine Laghaout
Date: 2014-12-17

Generate a performance report for recursions going from 1 to 'maxN' and
'etaSamples'+1 equally-spaced values of the quantum efficiency from 
'etaMin' to 'etaMax'. For all other input parameters, refer to main.m.

The performance report is saved in the local directory as a time-stamped 
*.mat file with the following (etaSamples+1)-by-maxN matrices:

D: Distinguishability of the potential states at the leaves.
avgOverlap: Average of the ratio of the minimum to maximum probabilities at
    each leaf
DiscErr: Probabiliy of error resulting of applying Pi1 and Pi2 to the
    input states
B: The Bell factor obtained from a tri-partite W state
mainClock: Timer of the whole function

Also recorded are the key input values.
%}
function Reporter(FockDim, Operation, ParamsBounds, Pi, M, C, minN, ...
    maxN, etaSamples, etaMin, etaMax, FoM_ID)

%ReporterClock = tic;
%ProgressMonitor = waitbar(0, 'Please wait...');
ReportID = ['Reporter', datestr(now, 30), '.mat'];

% Take into account the possibility that only one value of 'eta' is to be
% sampled and check that the array of sample quantum efficiencies is valid.
if etaSamples > 0 && etaMin < etaMax
    etaArray = etaMin:(etaMax-etaMin)/etaSamples:etaMax;
elseif etaSamples == 0 || etaMin == etaMax
    etaSamples = 0;
    etaArray = etaMin;
else
    error('Invalid array of samples for quantum efficiency.');
end

SEQ = NaN;
D = NaN(etaSamples+1, maxN);
avgOverlap = NaN(etaSamples+1, maxN);
DiscErr = NaN(etaSamples+1, maxN);
B = NaN(etaSamples+1, maxN);
mainClock = NaN(etaSamples+1, maxN);

save(ReportID, 'D', 'avgOverlap', 'DiscErr', 'B', 'mainClock', ...
    'Operation', 'ParamsBounds', 'Pi', 'M', 'C', 'minN', 'maxN', ...
    'etaMin', 'etaMax', 'etaSamples', 'FoM_ID');

for N = minN:maxN
    
    etaIndex = 1;

    for eta = etaArray

        [D(etaIndex, N), B(etaIndex, N), ~, ~, avgOverlap(etaIndex, N), ...
            DiscErr(etaIndex, N), mainClock(etaIndex, N)] = ...
            main(FockDim, Operation, ParamsBounds, Pi, eta, M, C, N, SEQ, FoM_ID);
        
        etaIndex = etaIndex + 1;
        PercentageDone = ((N-minN)*(etaSamples+1)+etaIndex-1)...
            /((etaSamples+1)*(maxN-minN+1));

        fprintf('%.2f%%\t', 100*PercentageDone);
        %fprintf('N = %d\teta = %.2f%%\n', N, 100*eta);
        %waitbar(PercentageDone, ProgressMonitor, sprintf('%d%% [%d s]', round(100*PercentageDone), round(toc(ReporterClock))));
       
    end

    save(ReportID, 'D', 'avgOverlap', 'DiscErr', 'B', 'mainClock', ...
        '-append');

end

fprintf('\n');

%close(ProgressMonitor);
