%% Conducts sensitivity analysis 

% Parameters tested: discount rate, development parameters, damage
% function, rate of sand cost increase, rate of storm induced erosion
% increase, robustness to SLR scenario, time horizon, non-structural value
% at risk

% Saves output to specified output folder

%% Specify location to save outputs
output_folder = 'OutputJun19_polynomial';

%% Specify tests to run
Discount = 0;
Development = 0; 
Damage = 0;
SandCost = 0;
Storms = 0;
SeaLevel = 0;
Horizon = 0;
W = 0;
tax = 0;

%% Tax Rate
taxVals = 5;
if tax
    pars = parameters(0); 
    tax_vec = linspace(pars.lMin,pars.lMax,taxVals);
    strategy = cell(taxVals,3);
    optS = cell(taxVals,3);
    NPV = cell(taxVals,3);
    for scenario = 0:2
        pars = parameters(scenario);
        for i = 1:taxVals
            pars.l = tax_vec(i);   
            
            % Solve system and save outputs of interest
            [S,actions,x,V,~,presentVal] = main(pars);
            optS(i,scenario+1) = {[x V S]};
            strategy(i,scenario+1) = {actions};
            NPV(i,scenario+1) = {presentVal};
        end
    end
    TaxOutput.NPV = NPV;
    TaxOutput.optS = optS;
    TaxOutput.strategy = strategy;
    TaxOutput.l = tax_vec;
    save(strcat(output_folder,'/TaxOutput'),'TaxOutput');
end


%% Other value
WVals = 5;
if W
    pars = parameters(0); 
    W_vec = linspace(pars.WMin,pars.WMax,WVals);
    strategy = cell(WVals,3);
    optS = cell(WVals,3);
    NPV = cell(WVals,3);
    for scenario = 0:2
        pars = parameters(scenario);
        for i = 1:WVals
            pars.W = W_vec(i);   
            
            % Solve system and save outputs of interest
            [S,actions,x,V,~,presentVal] = main(pars);
            optS(i,scenario+1) = {[x V S]};
            strategy(i,scenario+1) = {actions};
            NPV(i,scenario+1) = {presentVal};
        end
    end
    WOutput.NPV = NPV;
    WOutput.optS = optS;
    WOutput.strategy = strategy;
    WOutput.W = W_vec;
    save(strcat(output_folder,'/WOutput'),'WOutput');
end

%% Effect of the discount rate
if Discount
    pars = parameters(0); 
    deltas = [pars.deltaMin;pars.delta;pars.deltaMax];
    strategy = cell(length(deltas),3);
    optS = cell(length(deltas),3);
    NPV = cell(length(deltas),3);
    for scenario = 0:2
        pars = parameters(scenario);
        for i = 1:length(deltas)
            pars.delta = deltas(i);   
            
            % Solve system and save outputs of interest
            [S,actions,x,V,~,presentVal] = main(pars);
            optS(i,scenario+1) = {[x V S]};
            strategy(i,scenario+1) = {actions};
            NPV(i,scenario+1) = {presentVal};
        end
    end
    DiscountOutput.NPV = NPV;
    DiscountOutput.optS = optS;
    DiscountOutput.strategy = strategy;
    DiscountOutput.deltas = deltas;
    save(strcat(output_folder,'/DiscountOutput'),'DiscountOutput');
end

%% Effect of coastal development and property value parameters
if Development
    
    % Effect of development parameters: alpha and beta
    scenarios = 0:2; 
    NPV = cell(3,3,length(scenarios));
    optS = cell(3,3,length(scenarios));
    strategy = cell(3,3,length(scenarios));
    for s = 1:length(scenarios)
        scenario = scenarios(s);
        pars = parameters(scenario);
        alpha = linspace(pars.alphaMin,pars.alphaMax,21);
        beta = [pars.betaMin,pars.beta,pars.betaMax];
        for i = 1:length(alpha)
            pars.alpha = alpha(i);
            for j = 1:length(beta)
                pars.beta = beta(j);
                    [S,actions,x,V,~,presentVal] = main(pars);
                    optS(i,j,s) = {[x V S]};
                    strategy(i,j,s) = {actions};
                    NPV(i,j,s) = {presentVal};
            end
        end
    end
    DevelopmentOutput.NPV = NPV;
    DevelopmentOutput.optS = optS;
    DevelopmentOutput.strategy = strategy;
    DevelopmentOutput.alpha = alpha;
    DevelopmentOutput.beta = beta;
    save(strcat(output_folder,'/DevelopmentOutput'),'DevelopmentOutput');
    
    %Focus on alpha
    n1 = 15;
    NPV = cell(n1,3);
    optS = cell(n1,3);
    strategy = cell(n1,3);
    scenarios = [0 1 2];
    for s = 1:length(scenarios)
        scenario = scenarios(s);
        pars = parameters(scenario);
        alpha = linspace(1,1.35,n1);
        for i = 1:length(alpha)
            pars.alpha = alpha(i);
            [S,actions,x,V,~,presentVal] = main(pars);
            optS(i,s) = {[x V S]};
            strategy(i,s) = {actions};
            NPV(i,s) = {presentVal};
        end
    end
    AlphaOutput.NPV = NPV;
    AlphaOutput.optS = optS;
    AlphaOutput.strategy = strategy;
    AlphaOutput.alpha = alpha;
    save(strcat(output_folder,'/AlphaOutput'),'AlphaOutput');
    
    % Effect of exponential development rate
    n2 = 11;
    d = linspace(pars.dMin,pars.dMax,n2);
    strategy = cell(n2,3);
    optS = cell(n2,3);
    NPV = cell(n2,3);
    valueFunc = zeros(n2,3);
    for scenario = 0:2
        pars = parameters(scenario);
        for i = 1:n2
            pars.d = d(i);           
            % Solve system and save outputs of interest
            [S,actions,x,V,~,presentVal,~,~,v] = main(pars);
            optS(i,scenario+1) = {[x V S]};
            strategy(i,scenario+1) = {actions};
            NPV(i,scenario+1) = {presentVal};
            valueFunc(i,scenario+1) = v;
        end
    end
    dOutput.NPV = NPV;
    dOutput.optS = optS;
    dOutput.strategy = strategy;
    dOutput.d = d;
    dOutput.valueFunc = valueFunc;
    save(strcat(output_folder,'/dOutput'),'dOutput');
    
    % Dependence on baseline property value
    n3 = 15;
    A = linspace(pars.AMin,pars.AMax,n3);
    strategy = cell(n3,3);
    optS = cell(n3,3);
    NPV = cell(n3,3);
    valueFunc = zeros(n3,3);
    for scenario = 0:2
        pars = parameters(scenario);
        for i = 1:n3
            pars.A = A(i);           
            % Solve system and save outputs of interest
            [S,actions,x,V,~,presentVal,~,~,v] = main(pars);
            optS(i,scenario+1) = {[x V S]};
            strategy(i,scenario+1) = {actions};
            NPV(i,scenario+1) = {presentVal};
            valueFunc(i,scenario+1) = v;
        end
    end
    AOutput.NPV = NPV;
    AOutput.optS = optS;
    AOutput.strategy = strategy;
    AOutput.A = A;
    AOutput.valueFunc = valueFunc;
    save(strcat(output_folder,'/AOutput'),'AOutput');
end

%% Damage functions
if Damage
    funcs = {'polynomial','linear','exponential','concave'};
    n = length(funcs);
    NPV = cell(n,3);
    optS = cell(n,3);
    strategy = cell(n,3);
    valueFunc = zeros(n,3);
    scenarios = [0 1 2];
    for s = 1:length(scenarios)
        scenario = scenarios(s);
        pars = parameters(scenario); 
        for i = 1:n
            pars.Cfunc = funcs{i};
            [S,actions,x,V,~,presentVal,~,~,v] = main(pars);
            optS(i,s) = {[x V S]};
            strategy(i,s) = {actions};
            NPV(i,s) = {presentVal};
            valueFunc(i,s) = v;
        end
        
    end
    DamageOutput.NPV = NPV;
    DamageOutput.optS = optS;
    DamageOutput.strategy = strategy;
    DamageOutput.valueFunc = valueFunc;
    DamageOutput.functions = funcs;
    save(strcat(output_folder,'/DamageOutput'),'DamageOutput');
end

%% Sand cost
if SandCost
    n = 21;
    NPV = cell(n,3);
    optS = cell(n,3);
    strategy = cell(n,3);
    scenarios = [0 1 2];
    for s = 1:length(scenarios)
        scenario = scenarios(s);
        pars = parameters(scenario);
        xi = linspace(0,.1,n);
        for i = 1:n
            pars.xi = xi(i);
            [S,actions,x,V,~,presentVal] = main(pars);
            optS(i,s) = {[x V S]};
            strategy(i,s) = {actions};
            NPV(i,s) = {presentVal};
        end
    end
    SandCostOutput.NPV = NPV;
    SandCostOutput.optS = optS;
    SandCostOutput.strategy = strategy;
    SandCostOutput.xi = xi;
    save(strcat(output_folder,'/SandCostOutput'),'SandCostOutput');
end

%% Storm induced erosion
if Storms
    n = 21;
    NPV = cell(n,3);
    optS = cell(n,3);
    strategy = cell(n,3);
    scenarios = [0 1 2];
    for s = 1:length(scenarios)
        scenario = scenarios(s);
        pars = parameters(scenario);
        epsilon = linspace(0,1,n);
        for i = 1:n
            pars.epsilon = epsilon(i);
            [S,actions,x,V,~,presentVal] = main(pars);
            optS(i,s) = {[x V S]};
            strategy(i,s) = {actions};
            NPV(i,s) = {presentVal};
        end
    end
    StormOutput.NPV = NPV;
    StormOutput.optS = optS;
    StormOutput.strategy = strategy;
    StormOutput.epsilon = epsilon;
    save(strcat(output_folder,'/StormOutput'),'StormOutput');
end

%% Robustness to sea level rise
if SeaLevel

    pars0 = parameters(0);
    pars1 = parameters(1);
    pars2 = parameters(2); 
    slr_pars = {pars0;pars1;pars2};
    [~,actions0] = main(pars0);
    [~,actions1] = main(pars1);
    [~,actions2] = main(pars2);
    slr_actions = [actions0,actions1,actions2];
    NPV = cell(3,3);
    for i = 1:3
        for j = 1:3
            [~,~,presentVal] = coupledSystem(slr_pars{i},slr_actions(:,j));
            NPV(i,j) = {presentVal};
        end
    end
    SLROutput.NPV = NPV;
    SLROutput.strategies = slr_actions;
    SLROutput.parameters = slr_pars;
    save(strcat(output_folder,'/SLROutput'),'SLROutput');
end

%% Time horizon
if Horizon
    
    T = [10;20;50;80;100;Inf]; 
    NPV = cell(length(T),3);
    strategy = cell(length(T),3);
    optS = cell(length(T),3);
    for scenario = 0:2
        pars = parameters(scenario);
        for i = 1:length(T)
            pars.T = T(i);
            
            % Solve system and save outputs of interest
            [S,actions,x,V,~,presentVal] = main(pars);
            optS(i,scenario+1) = {[x V S]};
            strategy(i,scenario+1) = {actions};
            NPV(i,scenario+1) = {presentVal};
        end
    end
    HorizonOutput.NPV = NPV;
    HorizonOutput.optS = optS;
    HorizonOutput.strategy = strategy;
    HorizonOutput.timeHorizon = T;
    save(strcat(output_folder,'/HorizonOutput'),'HorizonOutput');
      
end
