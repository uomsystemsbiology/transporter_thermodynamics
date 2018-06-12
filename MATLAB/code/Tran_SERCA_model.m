% This function calculates the steady-state cycling rate for the Tran et
% al. (2009) SERCA model. The input to this function is SR calcium, with
% the other concentrations set to those in Figure 13 of Tran et al. (2009).
% The cycling rates are calculated for the simplified three-state model.

function v_ss = Tran_SERCA_model(Casr)
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;
    INIT_STATES = Casr;

    ALGEBRAIC = zeros(1,6);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, INIT_STATES);
    v_ss = ALGEBRAIC(1,6);
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component environment (second)');
    LEGEND_CONSTANTS(:,1) = strpad('k_p1 in component SERCA (second_order_rate_constant)');
    LEGEND_CONSTANTS(:,2) = strpad('k_p2 in component SERCA (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,3) = strpad('k_p3 in component SERCA (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,4) = strpad('k_m1 in component SERCA (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,5) = strpad('k_m2 in component SERCA (second_order_rate_constant)');
    LEGEND_CONSTANTS(:,6) = strpad('k_m3 in component SERCA (second_order_rate_constant)');
    LEGEND_CONSTANTS(:,7) = strpad('kdcai in component SERCA (millimolar)');
    LEGEND_CONSTANTS(:,8) = strpad('kdcasr in component SERCA (millimolar)');
    LEGEND_CONSTANTS(:,9) = strpad('kdh1 in component SERCA (millimolar)');
    LEGEND_CONSTANTS(:,10) = strpad('kdhi in component SERCA (millimolar_squared)');
    LEGEND_CONSTANTS(:,11) = strpad('kdhsr in component SERCA (millimolar_squared)');
    LEGEND_CONSTANTS(:,12) = strpad('kdh in component SERCA (millimolar)');
    LEGEND_CONSTANTS(:,13) = strpad('n in component SERCA (dimensionless)');
    LEGEND_CONSTANTS(:,14) = strpad('Ca_i in component SERCA (millimolar)');
    LEGEND_STATES(:,1) = strpad('Ca_sr in component SERCA (millimolar)');
    LEGEND_CONSTANTS(:,15) = strpad('H_i in component SERCA (millimolar)');
    LEGEND_CONSTANTS(:,16) = strpad('ATP in component SERCA (millimolar)');
    LEGEND_CONSTANTS(:,17) = strpad('ADP in component SERCA (millimolar)');
    LEGEND_CONSTANTS(:,18) = strpad('P_i in component SERCA (millimolar)');
    LEGEND_CONSTANTS(:,19) = strpad('T_Cai in component SERCA (dimensionless)');
    LEGEND_ALGEBRAIC(:,1) = strpad('T_Casr in component SERCA (dimensionless)');
    LEGEND_CONSTANTS(:,20) = strpad('T_H1 in component SERCA (dimensionless)');
    LEGEND_CONSTANTS(:,21) = strpad('T_Hi in component SERCA (dimensionless)');
    LEGEND_CONSTANTS(:,22) = strpad('T_Hsr in component SERCA (dimensionless)');
    LEGEND_CONSTANTS(:,23) = strpad('T_H in component SERCA (dimensionless)');
    LEGEND_CONSTANTS(:,24) = strpad('a_p1 in component SERCA (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,25) = strpad('a_p2 in component SERCA (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,2) = strpad('a_p3 in component SERCA (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,26) = strpad('a_m1 in component SERCA (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,3) = strpad('a_m2 in component SERCA (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,27) = strpad('a_m3 in component SERCA (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,4) = strpad('s1 in component SERCA (per_second_squared)');
    LEGEND_ALGEBRAIC(:,5) = strpad('s2 in component SERCA (per_second_squared)');
    LEGEND_CONSTANTS(:,28) = strpad('s3 in component SERCA (per_second_squared)');
    LEGEND_ALGEBRAIC(:,6) = strpad('v_cycle in component SERCA (first_order_rate_constant)');
    LEGEND_RATES(:,1) = strpad('d/dt Ca_sr in component SERCA (millimolar)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    CONSTANTS(:,1) = 25900;
    CONSTANTS(:,2) = 2540;
    CONSTANTS(:,3) = 20.5;
    CONSTANTS(:,4) = 2;
    CONSTANTS(:,5) = 67200;
    CONSTANTS(:,6) = 149;
    CONSTANTS(:,7) = 0.9;
    CONSTANTS(:,8) = 2.24;
    CONSTANTS(:,9) = 1.09e-5;
    CONSTANTS(:,10) = 3.54e-3;
    CONSTANTS(:,11) = 1.05e-8;
    CONSTANTS(:,12) = 7.24e-5;
    CONSTANTS(:,13) = 2;
    CONSTANTS(:,14) = 150e-6;
    STATES(:,1) = 0;
    CONSTANTS(:,15) = 1e-4;
    CONSTANTS(:,16) = 0.1;
    CONSTANTS(:,17) = 36.3e-3;
    CONSTANTS(:,18) = 15;
    CONSTANTS(:,19) = CONSTANTS(:,14)./CONSTANTS(:,7);
    CONSTANTS(:,28) = 1.00000;
    CONSTANTS(:,20) = CONSTANTS(:,15)./CONSTANTS(:,9);
    CONSTANTS(:,21) = power(CONSTANTS(:,15), CONSTANTS(:,13))./CONSTANTS(:,10);
    CONSTANTS(:,22) = power(CONSTANTS(:,15), CONSTANTS(:,13))./CONSTANTS(:,11);
    CONSTANTS(:,23) = CONSTANTS(:,15)./CONSTANTS(:,12);
    CONSTANTS(:,24) =  CONSTANTS(:,1).*CONSTANTS(:,16);
    CONSTANTS(:,25) = ( CONSTANTS(:,2).*power(CONSTANTS(:,19), 2.00000))./(power(CONSTANTS(:,19), 2.00000)+ power(CONSTANTS(:,19), 2.00000).*CONSTANTS(:,21)+ CONSTANTS(:,21).*(1.00000+CONSTANTS(:,20)));
    CONSTANTS(:,26) = ( CONSTANTS(:,4).*CONSTANTS(:,21))./(power(CONSTANTS(:,19), 2.00000)+ power(CONSTANTS(:,19), 2.00000).*CONSTANTS(:,21)+ CONSTANTS(:,21).*(1.00000+CONSTANTS(:,20)));
    CONSTANTS(:,27) =  CONSTANTS(:,6).*CONSTANTS(:,18);
    CONSTANTS(:,28) =  CONSTANTS(:,24).*CONSTANTS(:,25)+ CONSTANTS(:,27).*CONSTANTS(:,26)+ CONSTANTS(:,27).*CONSTANTS(:,25);
    if (isempty(STATES)), warning('Initial values for states not set');, end
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
    ALGEBRAIC(:,1) = STATES(:,1)./CONSTANTS(:,8);
    ALGEBRAIC(:,2) = ( CONSTANTS(:,3).*CONSTANTS(:,22))./( power(ALGEBRAIC(:,1), 2.00000).*CONSTANTS(:,23)+CONSTANTS(:,23)+ CONSTANTS(:,22).*(1.00000+CONSTANTS(:,23)));
    ALGEBRAIC(:,3) = ( CONSTANTS(:,5).*CONSTANTS(:,17).*power(ALGEBRAIC(:,1), 2.00000).*CONSTANTS(:,23))./( power(ALGEBRAIC(:,1), 2.00000).*CONSTANTS(:,23)+CONSTANTS(:,23)+ CONSTANTS(:,22).*(1.00000+CONSTANTS(:,23)));
    ALGEBRAIC(:,4) =  CONSTANTS(:,25).*ALGEBRAIC(:,2)+ CONSTANTS(:,26).*ALGEBRAIC(:,2)+ CONSTANTS(:,26).*ALGEBRAIC(:,3);
    ALGEBRAIC(:,5) =  CONSTANTS(:,24).*ALGEBRAIC(:,2)+ ALGEBRAIC(:,3).*CONSTANTS(:,24)+ ALGEBRAIC(:,3).*CONSTANTS(:,27);
    ALGEBRAIC(:,6) = ( CONSTANTS(:,24).*CONSTANTS(:,25).*ALGEBRAIC(:,2) -  CONSTANTS(:,26).*ALGEBRAIC(:,3).*CONSTANTS(:,27))./(ALGEBRAIC(:,4)+ALGEBRAIC(:,5)+CONSTANTS(:,28));
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

