% Creates struct containing all parameters, including minimum and maximum
% values needed for coastal system problems
%
% Input:
%   slr_scenario: sea level rise scenario. 
%       0 for low, 1 for intermediate, 2 for high
%
% Output:
%   pars: parameter struct object

function pars = parameters(slr_scenario)

% Sea level rise parameters
if slr_scenario~=0 && slr_scenario~=1 && slr_scenario~=2 && slr_scenario~=-1
    error('SLR scenario incorrectly specified')
end
pars.scenario = slr_scenario;
pars.a = 2.355e-3; % historical sea level change (m/yr)
if slr_scenario==0
    pars.b = 0; % low sea level rise accleration
elseif slr_scenario==1
    pars.b = 0.0271e-3; % intermediate sea level rise accleration
elseif slr_scenario==2
    pars.b = 0.113e-3; % high sea level rise accleration
else
    pars.b = 0;
    pars.a = 0;
end

% Beach and nourishment parameters
pars.x_nourished = 6.096; % nourished beach width (m)
pars.x_crit = 0; % beach width nourishment trigger (m)
pars.mu = 1;%(pars.x_nourished-pars.x_crit)/pars.x_nourished; % nourished portion of beach 
pars.init_rate = -0.0792; % historical shoreline change rate (m/yr)
pars.theta = 0.1; % exponential erosion rate
pars.r = 70.04; % slope of the active profile
pars.H = pars.init_rate+pars.r*pars.a; % Bruun Rule correction term (m)

% Initial conditions
pars.tau_init = 0; %initial years since nourishment
pars.v_init = 669000; % initial value at risk ($1000)
pars.x_init = pars.x_nourished; % initial beach width (m)

% Time parameters
pars.deltaT = 1; % time step (yr)
pars.Ti = 2020; % year of simulation start
pars.sim_length = 100; % simulation length (yr)
pars.Tswitch = 1; % Time when relocation becomes an option
pars.T = Inf; % Time horizon

% Expected storm induced erosion
pars.lambda = 0.35; % storm frequency
pars.m = 1.68; % GEV location parameter
pars.sigma = 4.24; % GEV scale parameter
pars.k = 0.277; % GEV shape parameter
meanGEV = pars.m+pars.sigma*(gamma(1-pars.k)-1)/pars.k;
p = poisspdf(1:4,pars.lambda);
pars.E = 0;
for n = 1:4
    M = 0;
    for i = 1:n
        M = M+meanGEV/i;
    end
    pars.E = pars.E + 0.1*p(n)*M; % annual expected storm induced erosion
end
pars.epsilon = 0; % increase in storm induced erosion with SLR

% Property value parameters
pars.d = 0; % development rate constant (1/yr)
pars.alpha = (1+0.01)^3.2808; % property value increase due to 1 m increase in beach width
pars.beta = 0.5;%(1-0.147)^3.2808; % property value decrease due to 1 m increase in sea level
pars.A = 669000;
pars.v_init = pars.A*(pars.alpha^pars.x_init); % baseline property value
pars.W = 5e5; % non-structural value at risk

% Benefit and cost parameters
pars.delta = 0.0275; % Discount rate
pars.eta = 824; % land value ($1000/m), assumes $14/sq ft and 5470 m of beach length
pars.l = 0.22;%0.0041; % St. Lucie County general fund tax rate
pars.c1 = 12000; % fixed cost of nourishment ($1000), assumes $14 million per nourishment, c2=350
pars.c2 = 350; % variable cost of nourishment ($1000/m), assumes $9.55/m^3, 5470 m of beach length, and 224,000 m^3 per 6.096 m nourishment
pars.xi = 0; % exponential increase in c2 as time progresses (0 means cost is autonomous)
pars.constructionCosts = 0;

pars.Cfunc = 'concave';
pars.phi_exp =  5.6999; % sea level base for proportion damaged
pars.phi_lin = 61.3951;
pars.phi_conc = 193.8357;
pars.phi_poly = 3.7625;
pars.kappa = 1.2; % beach width base for proportion damaged
pars.D0 = 5.4e-3; % expected proportion damaged when width=0 and sea level=0

% Relocation parameters
pars.relocationDelay = 1; % Years after decision is made that relocation occurs
pars.rho = 1; % Proportion of property value spent to relocate

% Feasibility constraints
pars.minInterval = 4;

%% max and min values for uncertainty and sensitivity analysis
pars.x_nourishedMin = 0.8*pars.x_nourished;
pars.x_nourishedMax = 1.2*pars.x_nourished;
pars.HMin = -0.2;
pars.HMax = 0.2;
pars.thetaMin = 0.8*pars.theta;
pars.thetaMax = 1.2*pars.theta;
pars.rMin = 0.8*pars.r;
pars.rMax = 1.2*pars.r;
pars.EMin = 0.8*pars.E;
pars.EMax = 1.2*pars.E;
pars.dMin = -0.05;
pars.dMax = 0.05;
pars.alphaMin = 1;
pars.alphaMax = 1.2;
pars.betaMin = .1;
pars.betaMax = 1;
pars.c1Min = 0.8*pars.c1;
pars.c1Max = 1.2*pars.c1;
pars.c2Min = 0.8*pars.c2;
pars.c2Max = 1.2*pars.c2;
pars.xiMin = 0;
pars.xiMax = 0.05;
pars.etaMin = 0.8*pars.eta;
pars.etaMax = 1.2*pars.eta;
pars.kappaMin = 1+0.8*(pars.kappa-1);
pars.kappaMax = 1+1.2*(pars.kappa-1);
pars.D0Min = 0.8*pars.D0;
pars.D0Max = 1.2*pars.D0;
pars.deltaMin = 0.01; 
pars.deltaMax = 0.07;
pars.relocationDelayMin = 1;
pars.relocationDelayMax = 10;
pars.bMin = 0;
pars.bMax = 0.113e-3;
pars.epsilonMin = 0;
pars.epsilonMax = 0.1;
pars.AMin = 0.1*pars.A;
pars.AMax = 3*pars.A;
pars.muMin = 0.5;
pars.muMax = 1;
pars.WMin = 2e5;
pars.WMax = 1e6;
pars.lMin = 0.176;
pars.lMax = 0.264;

end
