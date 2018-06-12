% This function solves the differential equation for a transporter model
% that couples the transport of a substrate to a second process that
% supplies energy

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = coupled_reactions_ode(tspan,init_states)
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;
    if exist('init_states','var')
        INIT_STATES = init_states;
    end

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

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
    LEGEND_CONSTANTS(:,4) = strpad('K_A in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,5) = strpad('K_B in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,6) = strpad('K_E1 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,7) = strpad('K_E2 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,8) = strpad('K_E3 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,9) = strpad('K_E4 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,10) = strpad('K_Se in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,11) = strpad('K_Si in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,12) = strpad('kappa_1 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,13) = strpad('kappa_2 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,14) = strpad('kappa_3 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,15) = strpad('kappa_4 in component environment (dimensionless)');
    LEGEND_STATES(:,1) = strpad('x_A in component environment (dimensionless)');
    LEGEND_STATES(:,2) = strpad('x_B in component environment (dimensionless)');
    LEGEND_STATES(:,3) = strpad('x_E1 in component environment (dimensionless)');
    LEGEND_STATES(:,4) = strpad('x_E2 in component environment (dimensionless)');
    LEGEND_STATES(:,5) = strpad('x_E3 in component environment (dimensionless)');
    LEGEND_STATES(:,6) = strpad('x_E4 in component environment (dimensionless)');
    LEGEND_STATES(:,7) = strpad('x_Se in component environment (dimensionless)');
    LEGEND_STATES(:,8) = strpad('x_Si in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,6) = strpad('e_1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,17) = strpad('e_2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,33) = strpad('e_3 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,19) = strpad('e_4 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,16) = strpad('e_5 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,20) = strpad('e_6 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,8) = strpad('e_7 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,22) = strpad('e_8 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,23) = strpad('e_9 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,2) = strpad('e_10 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,41) = strpad('e_11 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,24) = strpad('e_12 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,10) = strpad('e_13 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,26) = strpad('e_14 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,18) = strpad('e_15 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,47) = strpad('e_16 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,28) = strpad('e_17 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,12) = strpad('e_18 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,30) = strpad('e_19 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,58) = strpad('e_20 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,31) = strpad('e_21 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,14) = strpad('e_22 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,32) = strpad('e_23 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,4) = strpad('e_24 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,61) = strpad('f_1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,39) = strpad('f_2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,36) = strpad('f_3 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,38) = strpad('f_4 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,40) = strpad('f_5 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,37) = strpad('f_6 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,51) = strpad('f_7 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,48) = strpad('f_8 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,49) = strpad('f_9 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,52) = strpad('f_10 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,46) = strpad('f_11 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,45) = strpad('f_12 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,69) = strpad('f_13 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,67) = strpad('f_14 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,59) = strpad('f_15 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,57) = strpad('f_16 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,56) = strpad('f_17 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,72) = strpad('f_18 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,70) = strpad('f_19 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,68) = strpad('f_20 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,71) = strpad('f_21 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,73) = strpad('f_22 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,60) = strpad('f_23 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,62) = strpad('f_24 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,1) = strpad('mu_A in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,54) = strpad('v_A in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,3) = strpad('mu_B in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,64) = strpad('v_B in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,5) = strpad('mu_E1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,65) = strpad('v_E1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,7) = strpad('mu_E2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,55) = strpad('v_E2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,9) = strpad('mu_E3 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,74) = strpad('v_E3 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,11) = strpad('mu_E4 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,75) = strpad('v_E4 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,13) = strpad('mu_Se in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,76) = strpad('v_Se in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,15) = strpad('mu_Si in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,42) = strpad('v_Si in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,34) = strpad('Af_r1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,21) = strpad('Ar_r1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,35) = strpad('v_r1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,43) = strpad('Af_r2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,25) = strpad('Ar_r2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,44) = strpad('v_r2 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,27) = strpad('Af_r3 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,63) = strpad('Ar_r3 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,66) = strpad('v_r3 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,29) = strpad('Af_r4 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,50) = strpad('Ar_r4 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,53) = strpad('v_r4 in component environment (dimensionless)');
    LEGEND_RATES(:,1) = strpad('d/dt x_A in component environment (dimensionless)');
    LEGEND_RATES(:,2) = strpad('d/dt x_B in component environment (dimensionless)');
    LEGEND_RATES(:,3) = strpad('d/dt x_E1 in component environment (dimensionless)');
    LEGEND_RATES(:,4) = strpad('d/dt x_E2 in component environment (dimensionless)');
    LEGEND_RATES(:,5) = strpad('d/dt x_E3 in component environment (dimensionless)');
    LEGEND_RATES(:,6) = strpad('d/dt x_E4 in component environment (dimensionless)');
    LEGEND_RATES(:,7) = strpad('d/dt x_Se in component environment (dimensionless)');
    LEGEND_RATES(:,8) = strpad('d/dt x_Si in component environment (dimensionless)');
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
    CONSTANTS(:,10) = 1;
    CONSTANTS(:,11) = 1;
    CONSTANTS(:,12) = 1;
    CONSTANTS(:,13) = 1;
    CONSTANTS(:,14) = 1;
    CONSTANTS(:,15) = 1;
    STATES(:,1) = 100;
    STATES(:,2) = 1;
    STATES(:,3) = 0.5;
    STATES(:,4) = 0.5;
    STATES(:,5) = 0.5;
    STATES(:,6) = 0.5;
    STATES(:,7) = 100;
    STATES(:,8) = 10;
    CONSTANTS(:,16) = 0.000000;
    CONSTANTS(:,17) = 0.000000;
    CONSTANTS(:,18) = 0.000000;
    CONSTANTS(:,19) = 0.000000;
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
    RATES(:,1) = CONSTANTS(:,16);
    RATES(:,2) = CONSTANTS(:,17);
    RATES(:,7) = CONSTANTS(:,18);
    RATES(:,8) = CONSTANTS(:,19);
    ALGEBRAIC(:,5) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,6).*STATES(:,3));
    ALGEBRAIC(:,6) = ALGEBRAIC(:,5);
    ALGEBRAIC(:,17) = ALGEBRAIC(:,6);
    ALGEBRAIC(:,15) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,11).*STATES(:,8));
    ALGEBRAIC(:,16) = ALGEBRAIC(:,15);
    ALGEBRAIC(:,19) = ALGEBRAIC(:,16);
    ALGEBRAIC(:,33) = ALGEBRAIC(:,17)+ALGEBRAIC(:,19);
    ALGEBRAIC(:,34) = ALGEBRAIC(:,33);
    ALGEBRAIC(:,7) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,7).*STATES(:,4));
    ALGEBRAIC(:,8) = ALGEBRAIC(:,7);
    ALGEBRAIC(:,20) = ALGEBRAIC(:,8);
    ALGEBRAIC(:,21) = ALGEBRAIC(:,20);
    ALGEBRAIC(:,35) =  CONSTANTS(:,12).*(exp(ALGEBRAIC(:,34)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,21)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,37) = ALGEBRAIC(:,35);
    ALGEBRAIC(:,22) = ALGEBRAIC(:,8);
    ALGEBRAIC(:,1) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,4).*STATES(:,1));
    ALGEBRAIC(:,2) = ALGEBRAIC(:,1);
    ALGEBRAIC(:,23) = ALGEBRAIC(:,2);
    ALGEBRAIC(:,41) = ALGEBRAIC(:,22)+ALGEBRAIC(:,23);
    ALGEBRAIC(:,43) = ALGEBRAIC(:,41);
    ALGEBRAIC(:,9) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,8).*STATES(:,5));
    ALGEBRAIC(:,10) = ALGEBRAIC(:,9);
    ALGEBRAIC(:,24) = ALGEBRAIC(:,10);
    ALGEBRAIC(:,25) = ALGEBRAIC(:,24);
    ALGEBRAIC(:,44) =  CONSTANTS(:,13).*(exp(ALGEBRAIC(:,43)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,25)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,46) = ALGEBRAIC(:,44);
    ALGEBRAIC(:,48) = ALGEBRAIC(:,46);
    ALGEBRAIC(:,51) = ALGEBRAIC(:,37) - ALGEBRAIC(:,48);
    ALGEBRAIC(:,55) = ALGEBRAIC(:,51);
    RATES(:,4) = ALGEBRAIC(:,55);
    ALGEBRAIC(:,36) = ALGEBRAIC(:,35);
    ALGEBRAIC(:,39) = ALGEBRAIC(:,36);
    ALGEBRAIC(:,11) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,9).*STATES(:,6));
    ALGEBRAIC(:,12) = ALGEBRAIC(:,11);
    ALGEBRAIC(:,28) = ALGEBRAIC(:,12);
    ALGEBRAIC(:,29) = ALGEBRAIC(:,28);
    ALGEBRAIC(:,18) = ALGEBRAIC(:,6);
    ALGEBRAIC(:,3) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,5).*STATES(:,2));
    ALGEBRAIC(:,4) = ALGEBRAIC(:,3);
    ALGEBRAIC(:,32) = ALGEBRAIC(:,4);
    ALGEBRAIC(:,47) = ALGEBRAIC(:,18)+ALGEBRAIC(:,32);
    ALGEBRAIC(:,50) = ALGEBRAIC(:,47);
    ALGEBRAIC(:,53) =  CONSTANTS(:,15).*(exp(ALGEBRAIC(:,29)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,50)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,57) = ALGEBRAIC(:,53);
    ALGEBRAIC(:,59) = ALGEBRAIC(:,57);
    ALGEBRAIC(:,61) =  - ALGEBRAIC(:,39)+ALGEBRAIC(:,59);
    ALGEBRAIC(:,65) = ALGEBRAIC(:,61);
    RATES(:,3) = ALGEBRAIC(:,65);
    ALGEBRAIC(:,45) = ALGEBRAIC(:,44);
    ALGEBRAIC(:,26) = ALGEBRAIC(:,10);
    ALGEBRAIC(:,27) = ALGEBRAIC(:,26);
    ALGEBRAIC(:,30) = ALGEBRAIC(:,12);
    ALGEBRAIC(:,13) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,10).*STATES(:,7));
    ALGEBRAIC(:,14) = ALGEBRAIC(:,13);
    ALGEBRAIC(:,31) = ALGEBRAIC(:,14);
    ALGEBRAIC(:,58) = ALGEBRAIC(:,30)+ALGEBRAIC(:,31);
    ALGEBRAIC(:,63) = ALGEBRAIC(:,58);
    ALGEBRAIC(:,66) =  CONSTANTS(:,14).*(exp(ALGEBRAIC(:,27)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,63)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,67) = ALGEBRAIC(:,66);
    ALGEBRAIC(:,69) = ALGEBRAIC(:,45) - ALGEBRAIC(:,67);
    ALGEBRAIC(:,74) = ALGEBRAIC(:,69);
    RATES(:,5) = ALGEBRAIC(:,74);
    ALGEBRAIC(:,56) = ALGEBRAIC(:,53);
    ALGEBRAIC(:,68) = ALGEBRAIC(:,66);
    ALGEBRAIC(:,70) = ALGEBRAIC(:,68);
    ALGEBRAIC(:,72) =  - ALGEBRAIC(:,56)+ALGEBRAIC(:,70);
    ALGEBRAIC(:,75) = ALGEBRAIC(:,72);
    RATES(:,6) = ALGEBRAIC(:,75);
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
    ALGEBRAIC(:,5) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,6).*STATES(:,3));
    ALGEBRAIC(:,6) = ALGEBRAIC(:,5);
    ALGEBRAIC(:,17) = ALGEBRAIC(:,6);
    ALGEBRAIC(:,15) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,11).*STATES(:,8));
    ALGEBRAIC(:,16) = ALGEBRAIC(:,15);
    ALGEBRAIC(:,19) = ALGEBRAIC(:,16);
    ALGEBRAIC(:,33) = ALGEBRAIC(:,17)+ALGEBRAIC(:,19);
    ALGEBRAIC(:,34) = ALGEBRAIC(:,33);
    ALGEBRAIC(:,7) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,7).*STATES(:,4));
    ALGEBRAIC(:,8) = ALGEBRAIC(:,7);
    ALGEBRAIC(:,20) = ALGEBRAIC(:,8);
    ALGEBRAIC(:,21) = ALGEBRAIC(:,20);
    ALGEBRAIC(:,35) =  CONSTANTS(:,12).*(exp(ALGEBRAIC(:,34)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,21)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,37) = ALGEBRAIC(:,35);
    ALGEBRAIC(:,22) = ALGEBRAIC(:,8);
    ALGEBRAIC(:,1) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,4).*STATES(:,1));
    ALGEBRAIC(:,2) = ALGEBRAIC(:,1);
    ALGEBRAIC(:,23) = ALGEBRAIC(:,2);
    ALGEBRAIC(:,41) = ALGEBRAIC(:,22)+ALGEBRAIC(:,23);
    ALGEBRAIC(:,43) = ALGEBRAIC(:,41);
    ALGEBRAIC(:,9) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,8).*STATES(:,5));
    ALGEBRAIC(:,10) = ALGEBRAIC(:,9);
    ALGEBRAIC(:,24) = ALGEBRAIC(:,10);
    ALGEBRAIC(:,25) = ALGEBRAIC(:,24);
    ALGEBRAIC(:,44) =  CONSTANTS(:,13).*(exp(ALGEBRAIC(:,43)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,25)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,46) = ALGEBRAIC(:,44);
    ALGEBRAIC(:,48) = ALGEBRAIC(:,46);
    ALGEBRAIC(:,51) = ALGEBRAIC(:,37) - ALGEBRAIC(:,48);
    ALGEBRAIC(:,55) = ALGEBRAIC(:,51);
    ALGEBRAIC(:,36) = ALGEBRAIC(:,35);
    ALGEBRAIC(:,39) = ALGEBRAIC(:,36);
    ALGEBRAIC(:,11) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,9).*STATES(:,6));
    ALGEBRAIC(:,12) = ALGEBRAIC(:,11);
    ALGEBRAIC(:,28) = ALGEBRAIC(:,12);
    ALGEBRAIC(:,29) = ALGEBRAIC(:,28);
    ALGEBRAIC(:,18) = ALGEBRAIC(:,6);
    ALGEBRAIC(:,3) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,5).*STATES(:,2));
    ALGEBRAIC(:,4) = ALGEBRAIC(:,3);
    ALGEBRAIC(:,32) = ALGEBRAIC(:,4);
    ALGEBRAIC(:,47) = ALGEBRAIC(:,18)+ALGEBRAIC(:,32);
    ALGEBRAIC(:,50) = ALGEBRAIC(:,47);
    ALGEBRAIC(:,53) =  CONSTANTS(:,15).*(exp(ALGEBRAIC(:,29)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,50)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,57) = ALGEBRAIC(:,53);
    ALGEBRAIC(:,59) = ALGEBRAIC(:,57);
    ALGEBRAIC(:,61) =  - ALGEBRAIC(:,39)+ALGEBRAIC(:,59);
    ALGEBRAIC(:,65) = ALGEBRAIC(:,61);
    ALGEBRAIC(:,45) = ALGEBRAIC(:,44);
    ALGEBRAIC(:,26) = ALGEBRAIC(:,10);
    ALGEBRAIC(:,27) = ALGEBRAIC(:,26);
    ALGEBRAIC(:,30) = ALGEBRAIC(:,12);
    ALGEBRAIC(:,13) =  CONSTANTS(:,1).*CONSTANTS(:,2).*log( CONSTANTS(:,10).*STATES(:,7));
    ALGEBRAIC(:,14) = ALGEBRAIC(:,13);
    ALGEBRAIC(:,31) = ALGEBRAIC(:,14);
    ALGEBRAIC(:,58) = ALGEBRAIC(:,30)+ALGEBRAIC(:,31);
    ALGEBRAIC(:,63) = ALGEBRAIC(:,58);
    ALGEBRAIC(:,66) =  CONSTANTS(:,14).*(exp(ALGEBRAIC(:,27)./( CONSTANTS(:,1).*CONSTANTS(:,2))) - exp(ALGEBRAIC(:,63)./( CONSTANTS(:,1).*CONSTANTS(:,2))));
    ALGEBRAIC(:,67) = ALGEBRAIC(:,66);
    ALGEBRAIC(:,69) = ALGEBRAIC(:,45) - ALGEBRAIC(:,67);
    ALGEBRAIC(:,74) = ALGEBRAIC(:,69);
    ALGEBRAIC(:,56) = ALGEBRAIC(:,53);
    ALGEBRAIC(:,68) = ALGEBRAIC(:,66);
    ALGEBRAIC(:,70) = ALGEBRAIC(:,68);
    ALGEBRAIC(:,72) =  - ALGEBRAIC(:,56)+ALGEBRAIC(:,70);
    ALGEBRAIC(:,75) = ALGEBRAIC(:,72);
    ALGEBRAIC(:,38) = ALGEBRAIC(:,36);
    ALGEBRAIC(:,40) =  - ALGEBRAIC(:,38);
    ALGEBRAIC(:,42) = ALGEBRAIC(:,40);
    ALGEBRAIC(:,49) = ALGEBRAIC(:,46);
    ALGEBRAIC(:,52) =  - ALGEBRAIC(:,49);
    ALGEBRAIC(:,54) = ALGEBRAIC(:,52);
    ALGEBRAIC(:,60) = ALGEBRAIC(:,57);
    ALGEBRAIC(:,62) = ALGEBRAIC(:,60);
    ALGEBRAIC(:,64) = ALGEBRAIC(:,62);
    ALGEBRAIC(:,71) = ALGEBRAIC(:,68);
    ALGEBRAIC(:,73) = ALGEBRAIC(:,71);
    ALGEBRAIC(:,76) = ALGEBRAIC(:,73);
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

