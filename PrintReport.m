%{
Author: Amine Laghaout
Date: 2014-12-23

This function generates plots for the data saved by Reporter() and saves
them under the './figures/' directory.
%}
function PrintReport(matfile, name, N)

matfilePath = ['./reports/', matfile, '.mat'];

if exist(matfilePath, 'file') == 2

    MarkerSize = 8;
    FontSize = 23;
    load(matfilePath);
    
    if exist('N', 'var')
        minN = min(N);
        maxN = max(N);
    elseif ~exist('minN', 'var')
        minN = 1;
    end

    figure('visible', 'off');

    %% Bell factor
    
    if ~strcmp(func2str(Pi), 'HD')
        h_plot = plot(minN:maxN, B(1, :), ':xk', minN:maxN, B(2, :), '-.xk', minN:maxN, B(3, :), '-.dk', minN:maxN, B(4, :), '--sk', minN:maxN, B(5, :), '--ok', minN:maxN, B(6, :), '-ok', 'LineWidth', 3, 'MarkerSize', MarkerSize);
    else
        h_plot = plot(minN:maxN, B(1, :), '-ok', 'LineWidth', 3, 'MarkerSize', MarkerSize);        
    end
    xlabel('{\itN}', 'FontName', 'Times New Roman', 'FontSize', FontSize);
    %ylabel('Bell factor', 'FontName', 'Times New Roman', 'FontSize', FontSize);
    set(h_plot(:), 'LineWidth', 1.5);
    grid on;
    set(gca, 'ColorOrder', copper);
    xlhand = get(gca, 'xlabel');
    set(xlhand, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    ylhand = get(gca, 'ylabel');
    set(ylhand, 'FontName', 'Times New Roman', 'FontSize', FontSize);    
    axis([minN-0.2 maxN+0.2 min(min(B(:, minN:maxN)))-0.05 max(max(B(:, minN:maxN)))+0.05]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    if ~strcmp(func2str(Pi), 'HD') && false
        h_legend = legend('{\it\eta} = 50%', '{\it\eta} = 60%', '{\it\eta} = 70%', '{\it\eta} = 80%', '{\it\eta} = 90%', '{\it\eta} = 100%', 'Location', 'NorthOutside', 'Orientation', 'horizontal');
        set(h_legend, 'FontName', 'Times New Roman', 'FontSize', 11);    
    end        
    print('-dpng', ['.\figures\B_', name, '.png']);
    print('-depsc', ['.\figures\B_', name, '.eps']);

    %% Distinguishability

	if ~strcmp(func2str(Pi), 'HD')
        h_plot = plot(minN:maxN, 100*D(1, :), ':xk', minN:maxN, 100*D(2, :), '-.xk', minN:maxN, 100*D(3, :), '-.dk', minN:maxN, 100*D(4, :), '--sk', minN:maxN, 100*D(5, :), '--ok', minN:maxN, 100*D(6, :), '-ok', 'LineWidth', 3, 'MarkerSize', MarkerSize);
    else
        h_plot = plot(minN:maxN, 100*D(1, :), '-ok', 'LineWidth', 3, 'MarkerSize', MarkerSize);
    end
    xlabel('{\itN}', 'FontName', 'Times New Roman', 'FontSize', FontSize);
    %ylabel('Distinguishability [%]', 'FontName', 'Times New Roman', 'FontSize', FontSize);
    set(h_plot(:), 'LineWidth', 1.5);
    grid on;
    set(gca, 'ColorOrder', copper);
    xlhand = get(gca, 'xlabel');
    set(xlhand, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    ylhand = get(gca, 'ylabel');
    set(ylhand, 'FontName', 'Times New Roman', 'FontSize', FontSize);        
    axis([minN-0.2 maxN+0.2 100*max(0, min(min(D(:, minN:maxN)))-0.02) 100*(max(max(D(:, minN:maxN)))+0.02)]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    if ~strcmp(func2str(Pi), 'HD') && false
        h_legend = legend('{\it\eta} = 50%', '{\it\eta} = 60%', '{\it\eta} = 70%', '{\it\eta} = 80%', '{\it\eta} = 90%', '{\it\eta} = 100%', 'Location', 'NorthOutside', 'Orientation', 'horizontal');
        set(h_legend, 'FontName', 'Times New Roman', 'FontSize', 11);    
    end      
    print('-dpng', ['.\figures\D_', name, '.png']);    
    print('-depsc', ['.\figures\D_', name, '.eps']);    

    %% Discrimination error
    if ~strcmp(func2str(Pi), 'HD')
        h_plot = plot(minN:maxN, 100*DiscErr(1, :), ':xk', minN:maxN, 100*DiscErr(2, :), '-.xk', minN:maxN, 100*DiscErr(3, :), '-.dk', minN:maxN, 100*DiscErr(4, :), '--sk', minN:maxN, 100*DiscErr(5, :), '--ok', minN:maxN, 100*DiscErr(6, :), '-ok', 'LineWidth', 3, 'MarkerSize', MarkerSize);
    else
        h_plot = plot(minN:maxN, 100*DiscErr(1, :), '-ok', 'LineWidth', 3, 'MarkerSize', MarkerSize);   
    end    
    xlabel('{\itN}', 'FontName', 'Times New Roman', 'FontSize', FontSize);
    %ylabel('Discrimination error [%]', 'FontName', 'Times New Roman', 'FontSize', FontSize);
    set(h_plot(:), 'LineWidth', 1.5);
    grid on;
    set(gca, 'ColorOrder', copper);
    xlhand = get(gca, 'xlabel');
    set(xlhand, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    ylhand = get(gca, 'ylabel');
    set(ylhand, 'FontName', 'Times New Roman', 'FontSize', FontSize);        
    axis([minN-0.2 maxN+0.2 100*max(0, min(min(DiscErr(:, minN:maxN)))-0.02) 100*(max(max(DiscErr(:, minN:maxN)))+0.02)]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    if ~strcmp(func2str(Pi), 'HD') && false
        h_legend = legend('{\it\eta} = 50%', '{\it\eta} = 60%', '{\it\eta} = 70%', '{\it\eta} = 80%', '{\it\eta} = 90%', '{\it\eta} = 100%', 'Location', 'NorthOutside', 'Orientation', 'horizontal');
        set(h_legend, 'FontName', 'Times New Roman', 'FontSize', 11);    
    end     
    print('-dpng', ['.\figures\E_', name, '.png']);        
    print('-depsc', ['.\figures\E_', name, '.eps']);        

    %% Average overlap

    if ~strcmp(func2str(Pi), 'HD')
        h_plot = plot(minN:maxN, 100*avgOverlap(1, :), ':xk', minN:maxN, 100*avgOverlap(2, :), '-.xk', minN:maxN, 100*avgOverlap(3, :), '-.dk', minN:maxN, 100*avgOverlap(4, :), '--sk', minN:maxN, 100*avgOverlap(5, :), '--ok', minN:maxN, 100*avgOverlap(6, :), '-ok', 'LineWidth', 3, 'MarkerSize', MarkerSize);
    else
        h_plot = plot(minN:maxN, 100*avgOverlap(1, :), '-ok', 'LineWidth', 3, 'MarkerSize', MarkerSize); 
    end      
    xlabel('{\itN}', 'FontName', 'Times New Roman', 'FontSize', FontSize);
    %ylabel('Mean min-to-max ratio [%]', 'FontName', 'Times New Roman', 'FontSize', FontSize);
    set(h_plot(:), 'LineWidth', 1.5);
    grid on;
    set(gca, 'ColorOrder', copper);
    xlhand = get(gca, 'xlabel');
    set(xlhand, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    ylhand = get(gca, 'ylabel');
    set(ylhand, 'FontName', 'Times New Roman', 'FontSize', FontSize);        
    axis([minN-0.2 maxN+0.2 100*max(0, min(min(avgOverlap(:, minN:maxN)))-0.02) 100*(max(max(avgOverlap(:, minN:maxN)))+0.02)]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    if ~strcmp(func2str(Pi), 'HD') && false
        h_legend = legend('{\it\eta} = 50%', '{\it\eta} = 60%', '{\it\eta} = 70%', '{\it\eta} = 80%', '{\it\eta} = 90%', '{\it\eta} = 100%', 'Location', 'NorthOutside', 'Orientation', 'horizontal');
        set(h_legend, 'FontName', 'Times New Roman', 'FontSize', 11);    
    end
    print('-dpng', ['.\figures\R_', name, '.png']);      
    print('-depsc', ['.\figures\R_', name, '.eps']);      
    
    %{
    %% Runtime with its exponential fit

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
    h_plot = plot(minN:maxN, mainClock, '-ok', minN:maxN, fitClock, ':xk');
    xlabel('{\itN}', 'FontName', 'Times New Roman', 'FontSize', FontSize);
    %ylabel('Mean runtime [hrs]', 'FontName', 'Times New Roman', 'FontSize', FontSize);
    set(h_plot(:), 'LineWidth', 1.5);
    grid on;
    xlhand = get(gca, 'xlabel');
    set(xlhand, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    ylhand = get(gca, 'ylabel');
    set(ylhand, 'FontName', 'Times New Roman', 'FontSize', FontSize);        
    axis([minN-0.2 maxN+0.2 minClock maxClock]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    %%h_legend = legend('Actual', 'Fitted', 'Location', 'Best');
    %%set(h_legend, 'FontName', 'Times New Roman', 'FontSize', 11);  
    title(['b^{c\cdotN} with b = ', num2str(b), ', c = ', num2str(c)]);
    print('-dpng', ['.\figures\', matfile, '_mainClock.png']);

    %% Print summary
    
    fprintf('M = %d, C = %d, Operation = %s, Pi = %s\nTotal run time = %.2f hrs, ', ...
        M, C, func2str(Operation), func2str(Pi), mainClockTotal);
    fprintf('eta = [%.1f%%...%.1f%%], N = [%d...%d]\n', ...
        100*etaMin, 100*etaMax, minN, maxN);
    fprintf('ParamsBounds = [');
    fprintf('%f ', ParamsBounds);
    fprintf(']\n');    
    fprintf('Average runtime function b^(c*N): b = %f, c = %f\n', b, c);
    %}
else
    
    error('An invalid argument was passed.');
    
end