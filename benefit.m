% Calculates undiscounted benefits of beach and community location
%
% Inputs: 
%   S: nxm, matrix of states
%   pars: parameter struct with SLR and cost parameters
%
% Outputs:
%   B: total, undiscounted benefit n-vector

function B = benefit(S,pars)
[x,V] = xVL(S,pars);

B_beach = pars.eta*x; % recreational and habitat value of beach
B_location = pars.l*V; % value of having community in current location
B = B_beach + B_location; % total benefits ($1000)
end