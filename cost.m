% Calculates undiscounted costs of nourishment, relocation, and damages
%
% Inputs: 
%   S: matrix of states
%   A: strategy vector
%   pars: parameter struct with SLR and cost parameters
%
% Outputs:
%   C: total, undiscounted cost vector
%   nourishCost: total, undiscounted nourishment cost vector
%   relocateCost: total, undiscounted relocation cost vector
%   damageCost: total, undiscounted damage cost vector

function [C,nourishCost,relocateCost,damageCost] = cost(S,A,pars)

% Calculate beach width, property value, and sea level given current state
[x,V,L] = xVL(S,pars);
t = S(:,2);
R = S(:,3);
tau = S(:,1);

c2 = pars.c2*(1+pars.xi).^t;

nourishCost = (pars.c1+c2.*(pars.x_nourished-x)+pars.constructionCosts).*(A==1); % cost of nourishing the beach ($1000), make phi2 function of time?

relocateCost = pars.rho*V.*(R==pars.relocationDelay); % cost of relocating ($1000)

% Damage cost depends on damage function specified in pars
if strcmp(pars.Cfunc, 'linear')
    damageCost = (pars.D0*(1+L.*pars.phi_lin)./pars.kappa.^(x)).*(V+pars.W*(R<pars.relocationDelay));
elseif strcmp(pars.Cfunc, 'exponential')
    damageCost = (pars.D0*pars.phi_exp.^L./pars.kappa.^(x)).*(V+pars.W*(R<pars.relocationDelay));
elseif strcmp(pars.Cfunc, 'concave')
    damageCost = (pars.D0*(1+pars.phi_conc*(1-exp(-L)))./pars.kappa.^(x)).*(V+pars.W*(R<pars.relocationDelay));
else
    damageCost = (pars.D0*(1+L).^pars.phi_poly./pars.kappa.^(x)).*(V+pars.W*(R<pars.relocationDelay));
end

% Restrict nourishment interval to minimum of pars.minInterval
infeasible = (A==1) & (tau<(pars.minInterval-1));
feasibilityCost = zeros(size(infeasible));
feasibilityCost(infeasible) = Inf;

C = nourishCost + relocateCost + damageCost + feasibilityCost; % total cost ($1000)

end
