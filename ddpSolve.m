%% Solves discrete dynamic program for coastal systems
%
% Inputs:
%   pars: struct containing all parameter values
%   relocate: boolean to include relocation as option 
%       deafault = 1
%
% Outputs:
%   results: results struct from MDPSolve
%   model: model struct used to solve system
%
% Functions called:
%   cost
%   benefit
%   g
%
% Requires MDPSolve toolbox:
%   Fackler, P.L. 2017. MDPSOLVE Software for Dynamic Optimization.
%   Accessed 29 April 2018. https://sites.google.com/site/mdpsolve/home

function [results,model,S,X] = ddpSolve(pars,relocate)
%% Parameter values

S0 = [pars.tau_init,0,0,0];

if ~exist('relocate','var')
    relocate = 1;
end
if relocate
    A = [0,1,2]; %possible actions, 0:nothing, 1:nourish, 2:abandon
else
    A = [0,1]; %possible actions, 0:nothing, 1:nourish
end

%% State variables
tau = 0:pars.deltaT:pars.sim_length; %time since last nourishment, yr
time = 0:pars.deltaT:pars.sim_length;
relocation = 0:pars.deltaT:(pars.relocationDelay+1);
nourished = [0 1];

s = {tau',time',relocation',nourished'};
S = rectgrid(s); %state space

X = rectgrid(S,A'); %state-action grid
svars = 1:size(S,2);
Ix = getI(X,svars);

% Index of initial state
var0i = zeros(1,size(S,2));
for vari = 1:size(S,2)
    var0i(vari) = find(abs(s{vari}-S0(vari))==min(abs(s{vari}-S0(vari))),1,'last');
end
S0i = var0i(end);
place_factor = 1;
for vari = size(S,2)-1:-1:1
    place_factor = place_factor*length(s{vari+1});
    S0i = S0i + (var0i(vari)-1)*place_factor;
end

%% Reward Vector

costs = cost(X(:,svars),X(:,end),pars);
benefits = benefit(X(:,svars),pars);

reward = benefits-costs; %net benefits

%% Transition probability matrix
gX = @(X) g(X(:,svars),X(:,end),pars);
opt_g2p = struct('cleanup',2);
P = g2P(gX,s,X,[],[],opt_g2p);

%% Solve system
model.X = X;
model.Ix = Ix;
model.R = reward;
model.d = 1/(1+pars.delta);
model.P = P;
model.colstoch = 1;
model.svars = svars;
model.T = pars.T;

options.keepall = 1;

results = mdpsolve(model,options);

end
