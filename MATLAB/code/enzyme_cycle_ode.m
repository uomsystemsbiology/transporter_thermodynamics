% This function solves the differential equation for the enzyme cycle model
% of a transporter discussed in Figure 1.

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = enzyme_cycle_ode(tspan,init_states)
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;
    if exist('init_states','var')
        INIT_STATES = init_states;
    end

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 1);

    % Solve model with ODE solver
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI);
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('t in component environment (second)');
    LEGEND_CONSTANTS(:,1) = strpad('R in component environment (J_per_K_per_mol)');
    LEGEND_CONSTANTS(:,2) = strpad('T in component environment (kelvin)');
    LEGEND_CONSTANTS(:,3) = strpad('F in component environment (C_per_mol)');
    LEGEND_CONSTANTS(:,4) = strpad('K_E1 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,5) = strpad('K_E2 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,6) = strpad('K_Se in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,7) = strpad('K_Si in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,8) = strpad('kappa_1 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,9) = strpad('kappa_2 in component environment (dimensionless)');
    LEGEND_STATES(:,1) = strpad('x_E1 in component environment (dimensionless)');
    LEGEND_STATES(:,2) = strpad('x_E2 in component environment (dimensionless)');
    LEGEND_STATES(:,3) = strpad('x_Se in component environment (dimensionless)');
    LEGEND_STATES(:,4) = strpad('x_Si in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,17) = strpad('e_1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,9) = strpad('e_2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,2) = strpad('e_3 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,11) = strpad('e_4 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,15) = strpad('e_5 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,8) = strpad('e_6 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,4) = strpad('e_7 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,10) = strpad('e_8 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,16) = strpad('e_9 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,6) = strpad('e_10 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,25) = strpad('e_11 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,13) = strpad('e_12 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,20) = strpad('f_1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,23) = strpad('f_2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,34) = strpad('f_3 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,21) = strpad('f_4 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,22) = strpad('f_5 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,24) = strpad('f_6 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,32) = strpad('f_7 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,31) = strpad('f_8 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,33) = strpad('f_9 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,35) = strpad('f_10 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,29) = strpad('f_11 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,30) = strpad('f_12 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,1) = strpad('mu_E1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,37) = strpad('v_E1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,3) = strpad('mu_E2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,36) = strpad('v_E2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,5) = strpad('mu_Se in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,38) = strpad('v_Se in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,7) = strpad('mu_Si in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,26) = strpad('v_Si in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,18) = strpad('Af_r1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,12) = strpad('Ar_r1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,19) = strpad('v_r1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,14) = strpad('Af_r2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,27) = strpad('Ar_r2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,28) = strpad('v_r2 in component environment (dimensionless)');
    LEGEND_RATES(:,1) = strpad('d/dt x_E1 in component environment (dimensionless)');
    LEGEND_RATES(:,2) = strpad('d/dt x_E2 in component environment (dimensionless)');
    LEGEND_RATES(:,3) = strpad('d/dt x_Se in component environment (dimensionless)');
    LEGEND_RATES(:,4) = strpad('d/dt x_Si in component environment (dimensionless)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    CONSTANTS(:,1) = 8.314;
    CONSTANTS(:,2) = 310;
    CONSTANTS(:,3) = 96485;
    CONSTANTS(:,4) = 1;
    CONSTANTS(:,5) = 1;
    CONSTANTS(:,6) = 1;
    CONSTANTS(:,7) = 1;
    CONSTANTS(:,8) = 1;
    CONSTANTS(:,9) = 1;
    STATES(:,1) = 1;
    STATES(:,2) = 1;
    STATES(:,3) = 10;
    STATES(:,4) = 100;
    CONSTANTS(:,10) = 0.000000;
    CONSTANTS(:,11) = 0.000000;
    if (isempty(STATES)), warning('Initial values for states not set');, end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
        utilOnes = ones(statesRowCount, 1);
    end
    RATES(:,3) = CONSTANTS(:,10);
    RATES(:,4) = CONSTANTS(:,11);
    ALGEBRAIC(:,1) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,4).*STATES(:,1));
    ALGEBRAIC(:,2) = ALGEBRAIC(:,1);
    ALGEBRAIC(:,9) = ALGEBRAIC(:,2);
    ALGEBRAIC(:,7) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,7).*STATES(:,4));
    ALGEBRAIC(:,8) = ALGEBRAIC(:,7);
    ALGEBRAIC(:,15) = ALGEBRAIC(:,8);
    ALGEBRAIC(:,17) = ALGEBRAIC(:,9)+ALGEBRAIC(:,15);
    ALGEBRAIC(:,18) = ALGEBRAIC(:,17);
    ALGEBRAIC(:,3) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,5).*STATES(:,2));
    ALGEBRAIC(:,4) = ALGEBRAIC(:,3);
    ALGEBRAIC(:,11) = ALGEBRAIC(:,4);
    ALGEBRAIC(:,12) = ALGEBRAIC(:,11);
    ALGEBRAIC(:,19) =  CONSTANTS(:,8).*(exp(ALGEBRAIC(:,18)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,12)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,21) = ALGEBRAIC(:,19);
    ALGEBRAIC(:,13) = ALGEBRAIC(:,4);
    ALGEBRAIC(:,14) = ALGEBRAIC(:,13);
    ALGEBRAIC(:,10) = ALGEBRAIC(:,2);
    ALGEBRAIC(:,5) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,6).*STATES(:,3));
    ALGEBRAIC(:,6) = ALGEBRAIC(:,5);
    ALGEBRAIC(:,16) = ALGEBRAIC(:,6);
    ALGEBRAIC(:,25) = ALGEBRAIC(:,10)+ALGEBRAIC(:,16);
    ALGEBRAIC(:,27) = ALGEBRAIC(:,25);
    ALGEBRAIC(:,28) =  CONSTANTS(:,9).*(exp(ALGEBRAIC(:,14)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,27)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,30) = ALGEBRAIC(:,28);
    ALGEBRAIC(:,32) = ALGEBRAIC(:,21) - ALGEBRAIC(:,30);
    ALGEBRAIC(:,36) = ALGEBRAIC(:,32);
    RATES(:,2) = ALGEBRAIC(:,36);
    ALGEBRAIC(:,20) = ALGEBRAIC(:,19);
    ALGEBRAIC(:,23) = ALGEBRAIC(:,20);
    ALGEBRAIC(:,29) = ALGEBRAIC(:,28);
    ALGEBRAIC(:,31) = ALGEBRAIC(:,29);
    ALGEBRAIC(:,34) =  - ALGEBRAIC(:,23)+ALGEBRAIC(:,31);
    ALGEBRAIC(:,37) = ALGEBRAIC(:,34);
    RATES(:,1) = ALGEBRAIC(:,37);
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI)
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        utilOnes = ones(statesRowCount, 1);
    end
    ALGEBRAIC(:,1) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,4).*STATES(:,1));
    ALGEBRAIC(:,2) = ALGEBRAIC(:,1);
    ALGEBRAIC(:,9) = ALGEBRAIC(:,2);
    ALGEBRAIC(:,7) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,7).*STATES(:,4));
    ALGEBRAIC(:,8) = ALGEBRAIC(:,7);
    ALGEBRAIC(:,15) = ALGEBRAIC(:,8);
    ALGEBRAIC(:,17) = ALGEBRAIC(:,9)+ALGEBRAIC(:,15);
    ALGEBRAIC(:,18) = ALGEBRAIC(:,17);
    ALGEBRAIC(:,3) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,5).*STATES(:,2));
    ALGEBRAIC(:,4) = ALGEBRAIC(:,3);
    ALGEBRAIC(:,11) = ALGEBRAIC(:,4);
    ALGEBRAIC(:,12) = ALGEBRAIC(:,11);
    ALGEBRAIC(:,19) =  CONSTANTS(:,8).*(exp(ALGEBRAIC(:,18)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,12)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,21) = ALGEBRAIC(:,19);
    ALGEBRAIC(:,13) = ALGEBRAIC(:,4);
    ALGEBRAIC(:,14) = ALGEBRAIC(:,13);
    ALGEBRAIC(:,10) = ALGEBRAIC(:,2);
    ALGEBRAIC(:,5) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,6).*STATES(:,3));
    ALGEBRAIC(:,6) = ALGEBRAIC(:,5);
    ALGEBRAIC(:,16) = ALGEBRAIC(:,6);
    ALGEBRAIC(:,25) = ALGEBRAIC(:,10)+ALGEBRAIC(:,16);
    ALGEBRAIC(:,27) = ALGEBRAIC(:,25);
    ALGEBRAIC(:,28) =  CONSTANTS(:,9).*(exp(ALGEBRAIC(:,14)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,27)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,30) = ALGEBRAIC(:,28);
    ALGEBRAIC(:,32) = ALGEBRAIC(:,21) - ALGEBRAIC(:,30);
    ALGEBRAIC(:,36) = ALGEBRAIC(:,32);
    ALGEBRAIC(:,20) = ALGEBRAIC(:,19);
    ALGEBRAIC(:,23) = ALGEBRAIC(:,20);
    ALGEBRAIC(:,29) = ALGEBRAIC(:,28);
    ALGEBRAIC(:,31) = ALGEBRAIC(:,29);
    ALGEBRAIC(:,34) =  - ALGEBRAIC(:,23)+ALGEBRAIC(:,31);
    ALGEBRAIC(:,37) = ALGEBRAIC(:,34);
    ALGEBRAIC(:,22) = ALGEBRAIC(:,20);
    ALGEBRAIC(:,24) =  - ALGEBRAIC(:,22);
    ALGEBRAIC(:,26) = ALGEBRAIC(:,24);
    ALGEBRAIC(:,33) = ALGEBRAIC(:,29);
    ALGEBRAIC(:,35) = ALGEBRAIC(:,33);
    ALGEBRAIC(:,38) = ALGEBRAIC(:,35);
end

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end

