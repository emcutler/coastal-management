%% Calculates beach width, property value, and sea level given current state
%
% Inputs: 
%   S: vector of states
%   pars: parameter value struct
%   E: optional vector of storm induced erosion
%
% Outputs:
%   x: vector of beach widths
%   V: vector of property values
%   L: vector of sea levels

function [x,V,L] = xVL(S,pars,E)
tau = S(:,1);
t = S(:,2);
R = S(:,3);
nourishing = S(:,4);

L = pars.a*t+pars.b*(t.^2+2*t.*(pars.Ti-1992)); % sea level compared to beginning of simulation
L_int = 0.5*pars.a*t.^2 + pars.b*(t.^3/3+t.^2*(pars.Ti-1992)); %time integral of sea level rise since beginning of simulation

if ~exist('E','Var')
    E = pars.E*tau+pars.epsilon*L_int; % storm induced erosion
end

% Gradual erosion
gamma_erosion = zeros(size(tau)); 
for i = 1:length(tau)
    if tau(i)>0 && t(i)>0 && t(i)>=tau(i)
        L_t = pars.a*t(i)+pars.b*(t(i).^2+2*t(i).*(pars.Ti-1992));
        L_tau = pars.a*(t(i)-tau(i))+pars.b*((t(i)-tau(i)).^2+2*(t(i)-tau(i)).*(pars.Ti-1992));
        if tau(i)==t(i) % before first nourishment
            gamma_erosion(i) = pars.r*L_t - pars.H*t(i);
        else
            gamma_erosion(i) = pars.r*(L_t-L_tau) - pars.H*tau(i); %background erosion since last nourishment
        end
    else
        gamma_erosion(i) = 0;
    end
end

% Calculate x and V
x = pars.x_nourished.*nourishing + (max((1-pars.mu)*pars.x_nourished + pars.mu*exp(-pars.theta*tau)*pars.x_nourished - gamma_erosion - E,0)).*(1-nourishing);
V = (1+pars.d).^t.*pars.A.*(pars.alpha.^x).*(pars.beta.^L).*(R<=pars.relocationDelay);

end