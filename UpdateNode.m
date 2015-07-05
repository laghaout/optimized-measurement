%{
Author: Amine Laghaout
Date: 2014-12-15

A set of incoming states 'state' is mixed with the vacuum on a beam
splitter 'B'. The reflected mode undergoes transformation 'A' before
being projected on operator 'Pi'. The transmitted (i.e., updated) set of
states is returned as 'state'.
%}
function state = UpdateNode(state, B, A, Pi)

% Recover the dimension of the Fock Hilbert space
if length(A) == length(Pi) && length(B) == length(A)^2
    FockDim = length(A);
else
    error('Inconsistent matrix dimensions');
end

UFlow = 1e-12;
[sqrtPi, ~] = sqrtm(Pi);
Vacuum = zeros(FockDim);
Vacuum(1,1) = 1;

for i = 1:length(state) % For each potential state...

    % Incoming single mode state    
    rho = zeros(FockDim);
    rhoRange = 1:length(state(i).rho);
    rho(rhoRange, rhoRange) = state(i).rho;
    
    % Beam-splitting B, operation A, and projection Pi    
    rhotilde = kron(eye(FockDim), sqrtPi)*kron(eye(FockDim), A)*B*kron(rho, Vacuum)*ctranspose(B)*ctranspose(kron(eye(FockDim), A))*kron(eye(FockDim), sqrtPi);    
    
    % Trace out the reflected mode
    rho = PartialTrace(FockDim, rhotilde);
    
    P = trace(rho);
    
    % Make sure there is no underflow and check that the probability is 
    % otherwise valid
    if abs(P) < 1e-10 
        P = UFlow;
    elseif P < 0-UFlow || P > 1+UFlow || isnan(P) || abs(imag(P)) > 1e-6
        fprintf('Re(P) = %.21f, Im(P) = %.21f\n', real(P), imag(P));
        error('Invalid probability %s%% for state %d', 100*num2str(P), i);
    end    
    
    % Update the state and its classical probability
    state(i).rho = rho/P;
    state(i).P = state(i).P*P;
    
end
