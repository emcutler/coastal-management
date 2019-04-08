
%% Simulates shoreline change due to gradual erosion, storm impacts, and nourishment
%
% Inputs:
%   pars: struct object containing parameter values
%   strategy: vector containing information about actions
%       0: do nothing, 1: nourish, 2: relocate
%       default: nourish when x=x_crit, never relocate
%   stochastic: boolean indicating whether to use stochastic version.
%       default: 0 (deterministic). Set to 1 for stochastic system
%
% Outputs:
%   S: State matrix
%   strategy: strategy followed to generate S
%   NPV: total NPV following strategy
%   costs: vector of undiscounted costs
%   benefits: vector of undiscounted benefits
%   storms: struct containing number of storms and storm induced erosion
%
% Functions called:
%   cost
%   benefit

function [S,strategy,NPV,costs,benefits,storms] = coupledSystem(pars,strategy,stochastic)

%% Management strategy and sea level rise

% Default strategy if no strategy provided as input
if ~exist('strategy','Var')
    strategy = zeros(pars.sim_length,1);
    defaultStrat = 1;
else
    defaultStrat = 0;
end

S = zeros(pars.sim_length,6); 
relocateTime = find(strategy==2);

if ~isempty(relocateTime)
    S(relocateTime:relocateTime+pars.relocationDelay,3) = 1:pars.relocationDelay+1;
    S(relocateTime+1:end,3) = pars.relocationDelay+1;
end

S(:,2) = 0:pars.sim_length-1;
L = pars.a*S(:,2)+pars.b*(S(:,2).^2+2*S(:,2)*(pars.Ti-1992));

%% Poisson distribution for number of storms
p0 = poisscdf(0,pars.lambda);
p1 = poisscdf(1,pars.lambda);
p2 = poisscdf(2,pars.lambda);
p3 = poisscdf(3,pars.lambda);
p4 = poisscdf(4,pars.lambda);

%% Initialize variables
S(1,1) = pars.tau_init;

E = zeros(pars.sim_length,1); %storm induced shoreline change
n_storms = zeros(pars.sim_length,1); %number of storms

if ~exist('stochastic','Var')   
    stochastic = 0;
end
if stochastic
    %calculate number of storms and storm induced shoreline change
    rand_num = rand; 
    if rand_num < p0
        n_storms(1) = 0;
    elseif rand_num < p1
        n_storms(1) = 1;
    elseif rand_num < p2
        n_storms(1) = 2;
    elseif rand_num < p3
        n_storms(1) = 3;
    elseif rand_num < p4
        n_storms(1) = 4;
    else
        n_storms(1) = 5;
    end
    for i = 1:n_storms(1)
        E(1) = (E(1) + 0.1*gevrnd(pars.k,pars.sigma,pars.m)/i);
    end
else
    E(1) = pars.E*S(1,1) + pars.epsilon*L(1);
end



%% Run simulation
for t = 2:pars.sim_length
    if defaultStrat
        lastNourish = find(strategy==1,1,'last');
        if ~isempty(lastNourish)
            newE = sum(E(lastNourish:(t-1)));
        else
            newE = sum(E(1:(t-1)));
        end
        x = xVL(S(t-1,:),pars,newE);
        strategy(t-1) = x<=pars.x_crit;
    end
    if stochastic
        %Calculate number of storms in year t
        rand_num = rand; 
        if rand_num < p0
            n_storms(t) = 0;
        elseif rand_num < p1
            n_storms(t) = 1;
        elseif rand_num < p2
            n_storms(t) = 2;
        elseif rand_num < p3
            n_storms(t) = 3;
        elseif rand_num < p4
            n_storms(t) = 4;
        else
            n_storms(t) = 5;
        end        
        %Storm erosion using GEV distribution
        for i = 1:n_storms(t)
            E(t) = E(t) + 0.1*gevrnd(pars.k,pars.sigma,(pars.m+pars.epsilon*L(t-1)))/i;
        end
    else
        E(t) = pars.E + pars.epsilon*L(t-1);
    end
    S(t,1:4) = g(S(t-1,1:4),strategy(t-1),pars);    
    storms.n_storms = n_storms;
end

nourished_years = find(S(:,4)==1);

if ~isempty(nourished_years)
    E(1:(nourished_years(1)-1)) = cumsum(E(1:(nourished_years(1)-1)));
    for i = 1:length(nourished_years)-1
        E(nourished_years(i):(nourished_years(i+1)-1)) = cumsum(E(nourished_years(i):(nourished_years(i+1)-1)));
    end
    E(nourished_years(end):end) = cumsum(E(nourished_years(end):end));
else
    E = cumsum(E);
end
storms.E = E; % storm induced erosion

[x, V] = xVL(S(:,1:4),pars,E);
S(:,5:6) = [x, V];
%% Calculate net benefits
[C,nourishCost,relocateCost,damageCost] = cost(S(:,1:4),strategy,pars);
costs = [nourishCost,relocateCost,damageCost];
benefits = benefit(S(:,1:4),pars);
NPV = cumsum((benefits-C)./(1+pars.delta).^S(:,2));
end