% This function calculates the steady-state cycling rates for the
% Terkildsen et al. model, with updates to the equations and parameters as
% described in the manuscript. Membrane potential and extracellular Na are
% inputs to the function. The rest of the concentrations set to the values
% in Figure 9 of the manuscript

function v_ss = NaK_kinetic_ode(V,Nae)
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts(V,Nae);
    
    ALGEBRAIC = zeros(1,10);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, INIT_STATES);
    v_ss = ALGEBRAIC(1,10);

end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('t in component environment (second)');
    LEGEND_CONSTANTS(:,1) = strpad('R in component environment (mJ_per_K_per_mol)');
    LEGEND_CONSTANTS(:,2) = strpad('T in component environment (kelvin)');
    LEGEND_CONSTANTS(:,3) = strpad('F in component environment (C_per_mol)');
    LEGEND_CONSTANTS(:,4) = strpad('K0d_nak_nai in component environment (mM)');
    LEGEND_CONSTANTS(:,5) = strpad('K0d_nak_nae in component environment (mM)');
    LEGEND_CONSTANTS(:,6) = strpad('Kd_nak_nai2 in component environment (mM)');
    LEGEND_CONSTANTS(:,7) = strpad('Kd_nak_nae2 in component environment (mM)');
    LEGEND_CONSTANTS(:,8) = strpad('Kd_nak_ke in component environment (mM)');
    LEGEND_CONSTANTS(:,9) = strpad('Kd_nak_ki in component environment (mM)');
    LEGEND_CONSTANTS(:,10) = strpad('Kd_nak_MgATP in component environment (mM)');
    LEGEND_CONSTANTS(:,11) = strpad('rate_f_nak_1 in component environment (per_second)');
    LEGEND_CONSTANTS(:,12) = strpad('rate_f_nak_2 in component environment (per_second)');
    LEGEND_CONSTANTS(:,13) = strpad('rate_f_nak_3 in component environment (per_second)');
    LEGEND_CONSTANTS(:,14) = strpad('rate_f_nak_4 in component environment (per_second)');
    LEGEND_CONSTANTS(:,15) = strpad('rate_r_nak_1 in component environment (per_second_per_mM)');
    LEGEND_CONSTANTS(:,16) = strpad('rate_r_nak_2 in component environment (per_second)');
    LEGEND_CONSTANTS(:,17) = strpad('rate_r_nak_3 in component environment (per_second_per_mM2)');
    LEGEND_CONSTANTS(:,18) = strpad('rate_r_nak_4 in component environment (per_second)');
    LEGEND_CONSTANTS(:,19) = strpad('factor_nak_del in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,20) = strpad('pKd_nak_hpi in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,21) = strpad('Kd_nak_kpi in component environment (mM)');
    LEGEND_CONSTANTS(:,22) = strpad('Kd_nak_napi in component environment (mM)');
    LEGEND_STATES(:,1) = strpad('Vm in component environment (mV)');
    LEGEND_CONSTANTS(:,23) = strpad('Ki in component environment (mM)');
    LEGEND_CONSTANTS(:,24) = strpad('Ke in component environment (mM)');
    LEGEND_CONSTANTS(:,25) = strpad('Nai in component environment (mM)');
    LEGEND_CONSTANTS(:,26) = strpad('Nae in component environment (mM)');
    LEGEND_CONSTANTS(:,27) = strpad('MgATP in component environment (mM)');
    LEGEND_CONSTANTS(:,28) = strpad('MgADP in component environment (mM)');
    LEGEND_CONSTANTS(:,29) = strpad('Pi_total in component environment (mM)');
    LEGEND_CONSTANTS(:,30) = strpad('H in component environment (mM)');
    LEGEND_CONSTANTS(:,31) = strpad('Pi in component environment (mM)');
    LEGEND_ALGEBRAIC(:,1) = strpad('Kd_nak_nae1 in component environment (mM)');
    LEGEND_ALGEBRAIC(:,2) = strpad('Kd_nak_nai1 in component environment (mM)');
    LEGEND_ALGEBRAIC(:,3) = strpad('Nae_bound1 in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,4) = strpad('Nai_bound1 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,32) = strpad('Nae_bound2 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,33) = strpad('Nai_bound2 in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,34) = strpad('Ke_bound in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,35) = strpad('Ki_bound in component environment (dimensionless)');
    LEGEND_CONSTANTS(:,36) = strpad('MgATP_bound in component environment (dimensionless)');
    LEGEND_ALGEBRAIC(:,5) = strpad('alpha_f_nak_1 in component environment (per_second)');
    LEGEND_CONSTANTS(:,37) = strpad('alpha_f_nak_2 in component environment (per_second)');
    LEGEND_ALGEBRAIC(:,6) = strpad('alpha_f_nak_3 in component environment (per_second)');
    LEGEND_CONSTANTS(:,38) = strpad('alpha_f_nak_4 in component environment (per_second)');
    LEGEND_CONSTANTS(:,39) = strpad('alpha_r_nak_1 in component environment (per_second)');
    LEGEND_ALGEBRAIC(:,7) = strpad('alpha_r_nak_2 in component environment (per_second)');
    LEGEND_CONSTANTS(:,40) = strpad('alpha_r_nak_3 in component environment (per_second)');
    LEGEND_ALGEBRAIC(:,8) = strpad('alpha_r_nak_4 in component environment (per_second)');
    LEGEND_ALGEBRAIC(:,9) = strpad('sig_nak in component environment (per_second3)');
    LEGEND_ALGEBRAIC(:,10) = strpad('v_cyc_nak in component environment (per_second)');
    LEGEND_RATES(:,1) = strpad('d/dt Vm in component environment (mV)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts(V,Nae)
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    CONSTANTS(:,1) = 8314;
    CONSTANTS(:,2) = 310;
    CONSTANTS(:,3) = 96485;
    CONSTANTS(:,4) = 579.7295;
    CONSTANTS(:,5) = 0.0349;
    CONSTANTS(:,6) = 5.6399;
    CONSTANTS(:,7) = 1.0617e4;
    CONSTANTS(:,8) = 1.0817;
    CONSTANTS(:,9) = 1.6795e4;
    CONSTANTS(:,10) = 140.3709;
    CONSTANTS(:,11) = 1.4232e3;
    CONSTANTS(:,12) = 1.1565e4;
    CONSTANTS(:,13) = 194.4506;
    CONSTANTS(:,14) = 3.0630e4;
    CONSTANTS(:,15) = 225.9048;
    CONSTANTS(:,16) = 3.6355e4;
    CONSTANTS(:,17) = 2.8104e5;
    CONSTANTS(:,18) = 1.5740e6;
    CONSTANTS(:,19) = -0.0550;
    CONSTANTS(:,20) = 6.7700;
    CONSTANTS(:,21) = 292;
    CONSTANTS(:,22) = 224;
    STATES(:,1) = V;
    CONSTANTS(:,23) = 0;
    CONSTANTS(:,24) = 5.4;
    CONSTANTS(:,25) = 50;
    CONSTANTS(:,26) = Nae;
    CONSTANTS(:,27) = 10;
    CONSTANTS(:,28) = 0;
    CONSTANTS(:,29) = 0;
    CONSTANTS(:,30) = 3.9811e-5;
    CONSTANTS(:,31) = CONSTANTS(:,29)./(1.00000+CONSTANTS(:,23)./CONSTANTS(:,21)+CONSTANTS(:,30)./power(10.0000, 3.00000 - CONSTANTS(:,20))+CONSTANTS(:,25)./CONSTANTS(:,22));
    CONSTANTS(:,41) = 1.00000;
    CONSTANTS(:,32) = CONSTANTS(:,26)./CONSTANTS(:,7);
    CONSTANTS(:,33) = CONSTANTS(:,25)./CONSTANTS(:,6);
    CONSTANTS(:,34) = CONSTANTS(:,24)./CONSTANTS(:,8);
    CONSTANTS(:,35) = CONSTANTS(:,23)./CONSTANTS(:,9);
    CONSTANTS(:,36) = CONSTANTS(:,27)./CONSTANTS(:,10);
    CONSTANTS(:,37) = CONSTANTS(:,12);
    CONSTANTS(:,38) = ( CONSTANTS(:,14).*CONSTANTS(:,36))./(1.00000+CONSTANTS(:,36));
    CONSTANTS(:,39) =  CONSTANTS(:,15).*CONSTANTS(:,28);
    CONSTANTS(:,40) = ( CONSTANTS(:,17).*CONSTANTS(:,31).*CONSTANTS(:,30))./(1.00000+CONSTANTS(:,36));
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
    ALGEBRAIC(:,1) =  CONSTANTS(:,5).*exp(( (1.00000+CONSTANTS(:,19)).*CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2)));
    ALGEBRAIC(:,2) =  CONSTANTS(:,4).*exp(( CONSTANTS(:,19).*CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2)));
    ALGEBRAIC(:,3) = CONSTANTS(:,26)./ALGEBRAIC(:,1);
    ALGEBRAIC(:,4) = CONSTANTS(:,25)./ALGEBRAIC(:,2);
    ALGEBRAIC(:,5) = ( CONSTANTS(:,11).*ALGEBRAIC(:,4).*power(CONSTANTS(:,33), 2.00000))./(( ALGEBRAIC(:,4).*power(CONSTANTS(:,33), 2.00000)+power(1.00000+CONSTANTS(:,33), 2.00000)+power(1.00000+CONSTANTS(:,35), 2.00000)) - 1.00000);
    ALGEBRAIC(:,6) = ( CONSTANTS(:,13).*power(CONSTANTS(:,34), 2.00000))./(( ALGEBRAIC(:,3).*power(CONSTANTS(:,32), 2.00000)+power(1.00000+CONSTANTS(:,32), 2.00000)+power(1.00000+CONSTANTS(:,34), 2.00000)) - 1.00000);
    ALGEBRAIC(:,7) = ( CONSTANTS(:,16).*ALGEBRAIC(:,3).*power(CONSTANTS(:,32), 2.00000))./(( ALGEBRAIC(:,3).*power(CONSTANTS(:,32), 2.00000)+power(1.00000+CONSTANTS(:,32), 2.00000)+power(1.00000+CONSTANTS(:,34), 2.00000)) - 1.00000);
    ALGEBRAIC(:,8) = ( CONSTANTS(:,18).*power(CONSTANTS(:,35), 2.00000))./(( ALGEBRAIC(:,4).*power(CONSTANTS(:,33), 2.00000)+power(1.00000+CONSTANTS(:,33), 2.00000)+power(1.00000+CONSTANTS(:,35), 2.00000)) - 1.00000);
    ALGEBRAIC(:,9) =  CONSTANTS(:,39).*ALGEBRAIC(:,7).*CONSTANTS(:,40)+ CONSTANTS(:,39).*ALGEBRAIC(:,7).*CONSTANTS(:,38)+ CONSTANTS(:,39).*ALGEBRAIC(:,6).*CONSTANTS(:,38)+ CONSTANTS(:,37).*ALGEBRAIC(:,6).*CONSTANTS(:,38)+ ALGEBRAIC(:,7).*CONSTANTS(:,40).*ALGEBRAIC(:,8)+ ALGEBRAIC(:,5).*ALGEBRAIC(:,7).*CONSTANTS(:,40)+ ALGEBRAIC(:,5).*ALGEBRAIC(:,7).*CONSTANTS(:,38)+ ALGEBRAIC(:,5).*ALGEBRAIC(:,6).*CONSTANTS(:,38)+ CONSTANTS(:,39).*CONSTANTS(:,40).*ALGEBRAIC(:,8)+ CONSTANTS(:,37).*CONSTANTS(:,40).*ALGEBRAIC(:,8)+ ALGEBRAIC(:,5).*CONSTANTS(:,37).*CONSTANTS(:,40)+ ALGEBRAIC(:,5).*CONSTANTS(:,37).*CONSTANTS(:,38)+ CONSTANTS(:,39).*ALGEBRAIC(:,7).*ALGEBRAIC(:,8)+ CONSTANTS(:,39).*ALGEBRAIC(:,6).*ALGEBRAIC(:,8)+ CONSTANTS(:,37).*ALGEBRAIC(:,6).*ALGEBRAIC(:,8)+ ALGEBRAIC(:,5).*CONSTANTS(:,37).*ALGEBRAIC(:,6);
    ALGEBRAIC(:,10) = ( ALGEBRAIC(:,5).*CONSTANTS(:,37).*ALGEBRAIC(:,6).*CONSTANTS(:,38) -  CONSTANTS(:,39).*ALGEBRAIC(:,7).*CONSTANTS(:,40).*ALGEBRAIC(:,8))./ALGEBRAIC(:,9);
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

