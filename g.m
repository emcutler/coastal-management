%% State transition function
%
% Inputs: 
%   S: vector of states
%   A: vector of actions
%   pars: parameter value struct
%
% Outputs:
%   S: vector of states at next time period given S and A

function S = g(S,A,pars)

% Loop through all states in S and actions in A
for k = 1:size(S,1)
    S(k,2) = S(k,2) + pars.deltaT; %advance time
    
    if A(k) ~= 1 %not nourishing
        S(k,1) = S(k,1)+pars.deltaT; %advance tau
        S(k,4) = 0;
    else %nourishing
        S(k,1) = 0; %reset tau
        S(k,4) = 1;
    end
    
    if S(k,3) <= pars.relocationDelay && S(k,3) > 0 % planning to relocate 
        S(k,3) = S(k,3) + pars.deltaT; % advance relocation delay counter
    end

    if A(k) == 2 && S(k,3) == 0 % decide to relocate
        S(k,3) = pars.deltaT; % start relocation delay counter
    end
    
    S(k,1) = min(S(k,1),pars.sim_length);
    S(k,2) = min(S(k,2),pars.sim_length);
end
end