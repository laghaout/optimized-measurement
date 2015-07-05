%{
Author: Amine Laghaout
Date: 2015-02-06

INPUT

comvars: Cf. main.m
FockDim: Cf. main.m
M: Cf. main.m
state: Density matrix and classical probabilities of the potential states
B: Matrix of the beam splitter transformation in Fock basis
eta: Cf. main.m
Operation: Cf. main.m
Pi: Cf. main.m
ParamsBounds: Cf. main.m

OUTPUT

OptimalParam: Optimal parameter for the pre-projection operation
optFoM: Optimized figure of merit achieved locally at the current node
    (e.g. distinguishability, Bell factor)
Remarks: Remarks about the performance of the optimization procedure. For
    example, how many iterations were necessary? Was there convergence of
    the optimal results?
OptimalA: Operator of the optimal pre-projection operator corresponding to
    OptimalParam
%}
function [OptimalParam, optFoM, Remarks, OptimalA] = ...
    AdaptiveDecision(comvars, FockDim, M, state, B, ...
    eta, Operation, Pi, ParamsBounds, k)

% Recover the number of candidate states
C = length(state);  

% For each potential state we get a probability distribution function over
% the M possible outputs. It is these probability distribution functions
% between the different potential states that we are trying to get as
% distinguished as possible.
PDFs = NaN(C, M);

% String of possible remarks regarding the performance of the algorithm.
% E.g., the number of iterations, or whether the parameter bounds were
% reached without reaching an optimimum.
Remarks = '';

UFlow = 1e-15; % Threshold to prevent underflow errors

switch length(ParamsBounds)

    % The operation Â is to be determined at one predefined parameter only.
    case 1
        minParams = ParamsBounds;
        maxParams = ParamsBounds;

    % The operation Â is to be optimized over a range of numParameters 
    % that is specified by default.
    case 2
        minParams = ParamsBounds(1);
        maxParams = ParamsBounds(2);
        numParams = 10;

    % The operation Â is to be optimized over a range of numParameters 
    % that is passed as an argument.
    case 3
        minParams = ParamsBounds(1);
        maxParams = ParamsBounds(2);
        numParams = ParamsBounds(3);
        
    otherwise
        error('Invalid array of parameter bounds');
end

%% Find the optimal parameter between minParams and maxParams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if minParams < maxParams
    
    minConvergence = 1e-3; % Does does making this smaller worsen things?
    Iter = 5; % Set to 5 by default.
    OptimalParam = NaN;    
    OptimalA = NaN(FockDim);
    optFoMold = -Inf;
    optFoM = -1e100;
    
    % Normalize the local probabiliy distribution such that the
    % probabilities of each state at the node add up to unity. I.e. we
    % factor out any overall probability associated with the node.
    nodeP = sum([state(:).P]);
    for i = 1:C
        state(i).P = state(i).P/nodeP;
    end
    
    % There will be at most Iter iterations of the following loop unless it
    % is exited by a convergence smaller than minConvergence. At each
    % iteration of the loop, the array of parameters is refined and the
    % parity of the number of sample parameters is flipped.
    while abs(optFoM - optFoMold) > minConvergence && Iter > 0

        optFoMold = optFoM;
        Iter = Iter - 1;

        % Parameter sweep
        for Params = minParams:(maxParams-minParams)/numParams:maxParams

            A = Operation(FockDim, Params);

            % This switch block should be almost identically copied in the 
            % outer ie-elseif statment in the elseif minParams == maxParams 
            % case below. The only difference is that, there, we are only 
            % returning the figure of merit, we're not trying to see if 
            % it's being maximized.
            %%{
            switch comvars.FoM_ID
                
                %% Distinguishability-based optimization
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'D' 
                    
                    % Determine the local distinguishability
                    for m = 1:M
                        m_state = UpdateNode(state, B, A, Pi(FockDim, m, M, eta));
                        PDFs(:, m) = [m_state(:).P];
                    end
                    
                    FoM = Distinguishability(PDFs, [state(:).P]);

                    if FoM > optFoM
                        optFoM = FoM;
                        OptimalParam = Params;
                        OptimalA = A;
                    end 

                %% Bell factor-based optimization (Probability-independent)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'Bpi'

                    m1currKraus = comvars.Kraus; 
                    m2currKraus = comvars.Kraus; 
                    Pi1 = []; % Superoperator detecting the first qubit
                    Pi2 = []; % Superoperator detecting the second qubit

                    % WARNING: THIS IS HARD-CODED FOR TWO PROJECTORS ONLY
                    % AND THEREFORE TWO CANDIDATE STATES ONLY
                    
                    m = 1; % m = 1
                    [sqrtPi, ~] = sqrtm(Pi(FockDim, m, M, eta));
                    BlackBox = kron(eye(FockDim), sqrtPi)*kron(eye(FockDim), A)*B;               
                    for n = 0:FockDim-1
                        m1currKraus(k*FockDim+1:(k+1)*FockDim, n*FockDim+1:(n+1)*FockDim) = ...
                            ProjectMode2([n 0], 1, BlackBox);
                    end
                    
                    m = 2; % m = 2
                    [sqrtPi, ~] = sqrtm(Pi(FockDim, m, M, eta));
                    BlackBox = kron(eye(FockDim), sqrtPi)*kron(eye(FockDim), A)*B;               
                    for n = 0:FockDim-1
                        m2currKraus(k*FockDim+1:(k+1)*FockDim, n*FockDim+1:(n+1)*FockDim) = ...
                            ProjectMode2([n 0], 1, BlackBox);
                    end                    
                    
                    Pi1_Comb1 = cat(3, Pi1, m1currKraus);
                    Pi2_Comb1 = cat(3, Pi2, m2currKraus); 
                    Pi1_Comb2 = cat(3, Pi1, m2currKraus);
                    Pi2_Comb2 = cat(3, Pi2, m1currKraus);                     
                    
                    FoM_Comb1 = BellTest(FockDim, comvars.BellModes, ...
                        func2str(Pi), M, eta, Pi1_Comb1, Pi2_Comb1);
                    FoM_Comb2 = BellTest(FockDim, comvars.BellModes, ...
                        func2str(Pi), M, eta, Pi1_Comb2, Pi2_Comb2);
                    
                    FoM = max(FoM_Comb1, FoM_Comb2);

                    if FoM > optFoM
                        optFoM = FoM;
                        OptimalParam = Params;
                        OptimalA = A;
                    end
                    
                %% Bell factor-based optimization
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'B'

                    currKraus = comvars.Kraus; 
                    Pi1 = []; % Superoperator detecting the first qubit
                    Pi2 = []; % Superoperator detecting the second qubit

                    for m = 1:M % For each outcome...

                        % ... build the new Kraus operator...
                        
                        mPi = Pi(FockDim, m, M, eta);
                        m_state =  UpdateNode(state, B, A, mPi);
                        PDFs(:, m) = [m_state(:).P];
                        [sqrtPi, ~] = sqrtm(mPi);
                        BlackBox = kron(eye(FockDim), sqrtPi)*kron(eye(FockDim), A)*B;               
                        for n = 0:FockDim-1
                            currKraus(k*FockDim+1:(k+1)*FockDim, n*FockDim+1:(n+1)*FockDim) = ...
                                ProjectMode2([n 0], 1, BlackBox);
                        end

                        % ... and based on the probabilities of the
                        % candidate states, assign the new operator to
                        % either the `Pi1' or `Pi2' projector.
                        
                        [minProb, ~] = min(PDFs(:, m));
                        [maxProb, MostLikelyState] = max(PDFs(:, m));
                        if abs(maxProb-minProb) > UFlow
                            switch MostLikelyState
                                case 1
                                    Pi1 = cat(3, Pi1, currKraus);
                                case 2
                                    Pi2 = cat(3, Pi2, currKraus);
                                otherwise
                                    error('Invalid index of the most likely state at leaf %d', m);
                            end
                        else
                            Pi1 = cat(3, Pi1, currKraus);        
                        end      

                    end % Loop over M outcomes

                    FoM = BellTest(FockDim, comvars.BellModes, ...
                        func2str(Pi), M, eta, Pi1, Pi2);

                    if FoM > optFoM
                        optFoM = FoM;
                        OptimalParam = Params;
                        OptimalA = A;
                    end

                %% Discrimination error-based optimization
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'DiscErr'

                    % WARNING: WE ARE HARD-CODING THE SPECIAL CASE OF *TWO*
                    % EQUALLY SPACED AND EQUALLY WEIGHTED QUBITS.
                    rho1P = 1/2;
                    rho1 = zeros(FockDim);
                    rho1(1:2, 1:2) = Qubit(2*pi*(1-1)/C + pi/2, 0);
                    rho2P = 1/2;
                    rho2 = zeros(FockDim);
                    rho2(1:2, 1:2) = Qubit(2*pi*(2-1)/C + pi/2, 0); 

                    currKraus = comvars.Kraus; 
                    Pi1 = []; % Superoperator detecting the first qubit
                    Pi2 = []; % Superoperator detecting the second qubit

                    for m = 1:M % For each outcome...

                        % ... build the new Kraus operator...
                        
                        mPi = Pi(FockDim, m, M, eta);
                        m_state =  UpdateNode(state, B, A, mPi);
                        PDFs(:, m) = [m_state(:).P];
                        [sqrtPi, ~] = sqrtm(mPi);
                        BlackBox = kron(eye(FockDim), sqrtPi)*kron(eye(FockDim), A)*B;               
                        for n = 0:FockDim-1
                            currKraus(k*FockDim+1:(k+1)*FockDim, n*FockDim+1:(n+1)*FockDim) = ...
                                ProjectMode2([n 0], 1, BlackBox);
                        end                

                        % ... and based on the probabilities of the
                        % candidate states, assign the new operator to
                        % either the `Pi1' or `Pi2' projector.
                        
                        [minProb, ~] = min(PDFs(:, m));
                        [maxProb, MostLikelyState] = max(PDFs(:, m));
                        if abs(maxProb-minProb) > UFlow
                            switch MostLikelyState
                                case 1
                                    Pi1 = cat(3, Pi1, currKraus);
                                case 2
                                    Pi2 = cat(3, Pi2, currKraus);
                                otherwise
                                    error('Invalid index of the most likely state at leaf %d', m);
                            end
                        else
                            Pi1 = cat(3, Pi1, currKraus);        
                        end      

                    end % Loop over M outcomes                    

                    Pi1rho2 = trace(SuperOp(Pi1, rho2, 1, 1));
                    Pi2rho1 = trace(SuperOp(Pi2, rho1, 1, 1));
                    FoM = rho1P*Pi2rho1 + rho2P*Pi1rho2;            

                    if abs(FoM) < abs(optFoM)
                        optFoM = FoM;
                        OptimalParam = Params;
                        OptimalA = A;
                    end
                    
                %% Average min-to-max probability overlap optimization
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'avgOverlap'

                    % WARNING: THIS IS HARD-CODED FOR TWO CANDIDATE STATES
                    % ONLY
                    
                    FoM = 0;
                    
                    for m = 1:M
                        m_state = UpdateNode(state, B, A, Pi(FockDim, m, M, eta));
                        PDFs(:, m) = [m_state(:).P];
                        FoM = FoM + sum(PDFs(:, m))*min(PDFs(:, m))/max(PDFs(:, m));    
                    end                    

                    if abs(FoM) < abs(optFoM)
                        optFoM = FoM;
                        OptimalParam = Params;
                        OptimalA = A;
                    end                    
                    
                otherwise
                    
                    error('Invalid identifier for the figure of merit to be optimized.');
                    
            end % switch(comvars.FoM_ID)
            %}
            
        end % Inner parameter sweep (i.e., without convergence)

        % Re-centre the range of parameters if the optimal parameter is too
        % close to the end points of that range or if the range was preset
        % to be re-centered at OptimalParam via the flag
        % comvars.MOVABLE_RANGE. (FOR NOW, WE JUST IMPLEMENT THIS IF
        % comvars.MOVABLE_RANGE IS SET TO TRUE.)
        if comvars.MOVABLE_RANGE 
            %{
                || ...
                abs(maxParams - OptimalParam) < minConvergence || ...
                abs(OptimalParam - minParams) < minConvergence
            %}
                      
            if abs(OptimalParam - minParams) < minConvergence
                Remarks = 'Hit lower bound';
            elseif abs(maxParams - OptimalParam) < minConvergence
                Remarks = 'Hit upper bound';
            end

            minParams = OptimalParam - (maxParams-minParams)/2;
            maxParams = 2*OptimalParam - minParams;            
            
        end

        % Double the number of parameters and flip its parity so long as
        % there is more iterations to be carried out.
        if Iter > 0
            numParams = 2*numParams + 1 - mod(numParams, 2);
        end

    end % Outer parameter sweep (i.e., convergence loop)

    if ~isempty(Remarks)
        Remarks = [Remarks, '; '];
    end
    Remarks = [Remarks, num2str(numParams), ' samples'];

%% Evaluate at a fixed parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif minParams == maxParams
    
    OptimalParam = minParams;
    OptimalA = Operation(FockDim, OptimalParam);
    
    % This swicth statement should be almost identical to the one above.
    %%{
    switch comvars.FoM_ID

        %% Distinguishability-based optimization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'D' 

            % Determine the local distinguishability
            for m = 1:M
                m_state = UpdateNode(state, B, OptimalA, Pi(FockDim, m, M, eta));
                PDFs(:, m) = [m_state(:).P];
            end
            
            FoM = Distinguishability(PDFs, [state(:).P]);

        %% Bell factor-based optimization (Probability-independent)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'Bpi'

            m1currKraus = comvars.Kraus; 
            m2currKraus = comvars.Kraus; 
            Pi1 = []; % Superoperator detecting the first qubit
            Pi2 = []; % Superoperator detecting the second qubit

            % WARNING: THIS IS HARD-CODED FOR TWO PROJECTORS ONLY
            % AND THEREFORE TWO CANDIDATE STATES ONLY

            m = 1; % m = 1
            [sqrtPi, ~] = sqrtm(Pi(FockDim, m, M, eta));
            BlackBox = kron(eye(FockDim), sqrtPi)*kron(eye(FockDim), OptimalA)*B;               
            for n = 0:FockDim-1
                m1currKraus(k*FockDim+1:(k+1)*FockDim, n*FockDim+1:(n+1)*FockDim) = ...
                    ProjectMode2([n 0], 1, BlackBox);
            end

            m = 2; % m = 2
            [sqrtPi, ~] = sqrtm(Pi(FockDim, m, M, eta));
            BlackBox = kron(eye(FockDim), sqrtPi)*kron(eye(FockDim), OptimalA)*B;               
            for n = 0:FockDim-1
                m2currKraus(k*FockDim+1:(k+1)*FockDim, n*FockDim+1:(n+1)*FockDim) = ...
                    ProjectMode2([n 0], 1, BlackBox);
            end                    

            Pi1_Comb1 = cat(3, Pi1, m1currKraus);
            Pi2_Comb1 = cat(3, Pi2, m2currKraus); 
            Pi1_Comb2 = cat(3, Pi1, m2currKraus);
            Pi2_Comb2 = cat(3, Pi2, m1currKraus);                     

            FoM_Comb1 = BellTest(FockDim, comvars.BellModes, ...
                func2str(Pi), M, eta, Pi1_Comb1, Pi2_Comb1);
            FoM_Comb2 = BellTest(FockDim, comvars.BellModes, ...
                func2str(Pi), M, eta, Pi1_Comb2, Pi2_Comb2);

            FoM = max(FoM_Comb1, FoM_Comb2);                
                
        %% Bell factor-based optimization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'B'

            currKraus = comvars.Kraus; 
            Pi1 = []; % Superoperator detecting the first qubit
            Pi2 = []; % Superoperator detecting the second qubit

            for m = 1:M % For each outcome...

                % ... build the new Kraus operator...

                mPi = Pi(FockDim, m, M, eta);
                m_state =  UpdateNode(state, B, OptimalA, mPi);
                PDFs(:, m) = [m_state(:).P];
                [sqrtPi, ~] = sqrtm(mPi);
                BlackBox = kron(eye(FockDim), sqrtPi)*kron(eye(FockDim), OptimalA)*B;               
                for n = 0:FockDim-1
                    currKraus(k*FockDim+1:(k+1)*FockDim, n*FockDim+1:(n+1)*FockDim) = ...
                        ProjectMode2([n 0], 1, BlackBox);
                end                

                % ... and based on the probabilities of the
                % candidate states, assign the new operator to
                % either the `Pi1' or `Pi2' projector.

                [minProb, ~] = min(PDFs(:, m));
                [maxProb, MostLikelyState] = max(PDFs(:, m));
                if abs(maxProb-minProb) > UFlow
                    switch MostLikelyState
                        case 1
                            Pi1 = cat(3, Pi1, currKraus);
                        case 2
                            Pi2 = cat(3, Pi2, currKraus);
                        otherwise
                            error('Invalid index of the most likely state at leaf %d', m);
                    end
                else
                    Pi1 = cat(3, Pi1, currKraus);        
                end      

            end % Loop over M outcomes

            FoM = BellTest(FockDim, comvars.BellModes, ...
                func2str(Pi), M, eta, Pi1, Pi2);

        %% Discrimination error-based optimization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'DiscErr'

            % WARNING: WE ARE HARD-CODING THE SPECIAL CASE OF *TWO*
            % EQUALLY SPACED AND EQUALLY WEIGHTED QUBITS.
            rho1P = 1/2;
            rho1 = zeros(FockDim);
            rho1(1:2, 1:2) = Qubit(2*pi*(1-1)/C + pi/2, 0);
            rho2P = 1/2;
            rho2 = zeros(FockDim);
            rho2(1:2, 1:2) = Qubit(2*pi*(2-1)/C + pi/2, 0); 

            currKraus = comvars.Kraus; 
            Pi1 = []; % Superoperator detecting the first qubit
            Pi2 = []; % Superoperator detecting the second qubit

            for m = 1:M % For each outcome...

                % ... build the new Kraus operator...

                mPi = Pi(FockDim, m, M, eta);
                m_state =  UpdateNode(state, B, OptimalA, mPi);
                PDFs(:, m) = [m_state(:).P];
                [sqrtPi, ~] = sqrtm(mPi);
                BlackBox = kron(eye(FockDim), sqrtPi)*kron(eye(FockDim), OptimalA)*B;               
                for n = 0:FockDim-1
                    currKraus(k*FockDim+1:(k+1)*FockDim, n*FockDim+1:(n+1)*FockDim) = ...
                        ProjectMode2([n 0], 1, BlackBox);
                end                

                % ... and based on the probabilities of the
                % candidate states, assign the new operator to
                % either the `Pi1' or `Pi2' projector.

                [minProb, ~] = min(PDFs(:, m));
                [maxProb, MostLikelyState] = max(PDFs(:, m));
                if abs(maxProb-minProb) > UFlow
                    switch MostLikelyState
                        case 1
                            Pi1 = cat(3, Pi1, currKraus);
                        case 2
                            Pi2 = cat(3, Pi2, currKraus);
                        otherwise
                            error('Invalid index of the most likely state at leaf %d', m);
                    end
                else
                    Pi1 = cat(3, Pi1, currKraus);        
                end      

            end % Loop over M outcomes                    

            Pi1rho2 = trace(SuperOp(Pi1, rho2, 1, 1));
            Pi2rho1 = trace(SuperOp(Pi2, rho1, 1, 1));
            FoM = rho1P*Pi2rho1 + rho2P*Pi1rho2;            
            
        %% Average min-to-max probability overlap optimization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'avgOverlap'

            % WARNING: THIS IS HARD-CODED FOR TWO CANDIDATE STATES
            % ONLY

            FoM = 0;

            for m = 1:M
                m_state = UpdateNode(state, B, OptimalA, Pi(FockDim, m, M, eta));
                PDFs(:, m) = [m_state(:).P];
                FoM = FoM + sum(PDFs(:, m))*min(PDFs(:, m))/max(PDFs(:, m));    
            end                    

        otherwise

            error('Invalid identifier for the figure of merit to be optimized.');

    end % switch(comvars.FoM_ID)
    %}
    
    optFoM =  FoM;

%% Invalid parameter range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    
    error('Invalid parameter range');
    
end
