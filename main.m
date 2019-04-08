%% Main function used to solve system and simulate outputs
%
% Inputs:
%   pars: struct containing all parameter values
%       default: values set by parameters.m with low SLR
%
% Outputs:
%   optS: optimal states
%   actions: optimal strategy
%   NPV: net present value of optimal strategy
%   C: vector of undiscounted costs
%   benefits: vector of undiscounted benefits
%   valueFunc: optimal value function
%
% Functions called:
%   parameters (if pars unset)
%   ddpSolve
%   cost
%   benefit

function [optS,actions,NPV,C,B,valueFunc] = main(pars)

%% Parameter values
if ~exist('pars','var')
    pars = parameters(0);
    disp('default parameter values with low SLR used');
end
Tswitch = pars.Tswitch;

%% Solve dynamic program with and/or without relocation
if Tswitch > pars.sim_length %never include relocation
    [results1,~,S,X1] = ddpSolve(pars,0); %does not include relocation
elseif Tswitch == 1 %always include relocation
    [results2,~,S,X2] = ddpSolve(pars,1); %includes relocation
else %include relocation only after Tswitch
    [results1,~,S,X1] = ddpSolve(pars,0); %does not include relocation
    [results2,~,~,X2] = ddpSolve(pars,1); %includes relocation
end

%% Extract optimal states, actions, and NPV

S0i = find(S(:,1)==pars.tau_init,1); % index of initial state

actions = zeros(pars.sim_length,1);

if Tswitch > pars.sim_length
    valueFunc = results1.v(S0i);
    if isinf(pars.T)
        [Si,Xi] = mdpsim_inf(S0i,S,results1.Ixopt,X1,pars);
    else
        Si = S0i;
        Xi = [];
        for i = 1:ceil(pars.sim_length/pars.T)
            [nextSi,nextXi] = mdpsim_fin(Si(end),S,results1.Ixopt,X1,pars);
            if iscolumn(Si(1:end-1))
                Si = cat(1,Si(1:end-1),nextSi);
            else
                Si = cat(1,Si(1:end-1)',nextSi);
            end
            Xi = cat(1,Xi,nextXi);
        end
        Si = Si(1:pars.sim_length);
        Xi = Xi(1:pars.sim_length);
    end
    actions = X1(Xi,end);
elseif Tswitch == 1
    valueFunc = results2.v(S0i);
    if isinf(pars.T)
       [Si,Xi] = mdpsim_inf(S0i,S,results2.Ixopt,X2,pars);
    else
        Si = S0i;
        Xi = [];
        for i = 1:ceil(pars.sim_length/pars.T)
            [nextSi,nextXi] = mdpsim_fin(Si(end),S,results2.Ixopt,X2,pars);
            if iscolumn(Si(1:end-1))
                Si = cat(1,Si(1:end-1),nextSi);
            else
                Si = cat(1,Si(1:end-1)',nextSi);
            end
            Xi = cat(1,Xi,nextXi);
        end
        Si = Si(1:pars.sim_length);
        Xi = Xi(1:pars.sim_length);
    end
    actions = X2(Xi,end);
else
    valueFunc = NaN;
    if isinf(pars.T)
        [S1i,X1i] = mdpsim_inf(S0i,S,results1.Ixopt,X1,pars);
        Si(1:pars.Tswitch) = S1i(1:pars.Tswitch);
        
        [S2i,X2i] = mdpsim_inf(Si(pars.Tswitch),S,results2.Ixopt,X2,pars);
        Si(pars.Tswitch+1:pars.sim_length) = S2i(2:pars.sim_length-pars.Tswitch+1);
        
        actions(1:pars.Tswitch) = X1(X1i(1:pars.Tswitch),end);
        actions(pars.Tswitch+1:end) = X2(X2i(2:pars.sim_length-pars.Tswitch+1),end);
    else
        Si = S0i;
        X1i = [];
        for i = 1:ceil(pars.Tswitch/pars.T)
            [nextSi,nextXi] = mdpsim_fin(Si(end),S,results1.Ixopt,X1,pars);
            if iscolumn(Si(1:end-1))
                Si = cat(1,Si(1:end-1),nextSi);
            else
                Si = cat(1,Si(1:end-1)',nextSi);
            end
            X1i = cat(1,X1i,nextXi);
        end
        Si = Si(1:pars.Tswitch);
        X2i = [];
        for i = 1:ceil((pars.sim_length-pars.Tswitch)/pars.T)
            [nextSi,nextXi] = mdpsim_fin(Si(end),S,results2.Ixopt,X2,pars);
            if iscolumn(Si(1:end-1))
                Si = cat(1,Si(1:end-1),nextSi);
            else
                Si = cat(1,Si(1:end-1)',nextSi);
            end
            X2i = cat(1,X2i,nextXi);
        end
        Si = Si(1:pars.sim_length);
        actions(1:pars.Tswitch) = X1(X1i(1:pars.Tswitch),end);
        actions(pars.Tswitch+1:end) = X2(X2i(2:pars.sim_length-pars.Tswitch+1),end);
    end  
end

optS = S(Si,:);

%% Calculate costs and benefits
C = cost(optS,actions,pars);
B = benefit(optS,pars);
NPV = cumsum((B-C)./(1+pars.delta).^optS(:,2));

end

%%
% Simulate optimal time path for finite or infinite horizon problems

function [Si,Xi] = mdpsim_fin(S0i,S,Ix,X,pars)
T = pars.T;
Si = zeros(T+1,1);
Xi = zeros(T,1);
Si(1) = S0i;
Xi(1) = Ix(Si(1),1);
for t=2:T
    nextS = g(S(Si(t-1),1:end-1),X(Xi(t-1),end),pars);
    Si(t) = find(ismember(S,nextS,'rows'));
    Xi(t) = Ix(Si(t),t);
end
nextS = g(S(Si(T),1:end-1),X(Xi(T),end),pars);
Si(T+1) = find(ismember(S,nextS,'rows'));
end

function [Si,Xi] = mdpsim_inf(S0i,S,Ix,X,pars)
Si = zeros(pars.sim_length,1);
Xi = zeros(pars.sim_length,1);
Si(1) = S0i;
Xi(1) = Ix(Si(1));
for t = 2:pars.sim_length
    nextS = g(S(Si(t-1),1:end-1),X(Xi(t-1),end),pars);
    Si(t) = find(ismember(S,nextS,'rows'));
    Xi(t) = Ix(Si(t));
end
end