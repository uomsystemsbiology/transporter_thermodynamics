% This script simulates the bond graph of the Tran et al. (2009) SERCA
% model, and generates the plots for Figure 5.

clear;
clc;
close all;

%% Set directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep];
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];

%% Options
save_figures = true;

%% Define constants
R = 8.314; % Unit J/mol/K
T = 310; % Unit J/mol/K

W_i = 38.0; % Intracellular volume, unit pL
W_sr = 2.28; % SR volume, Unit pL
W_isr = W_i + W_sr; % Sum of intracellular and SR volumes, unit pL

%% Run system to steady state for different concentrations of Ca_sr
Casr_vec = transpose(0.05:0.05:5);
num_sim = length(Casr_vec);

v_ss_vec = zeros(num_sim,1);
v_ss_Tran = zeros(num_sim,1);

DG_vec = zeros(num_sim,1);

eff_vec = zeros(num_sim,1);

tspan = [0 5];

for i_sim = 1:num_sim
    % Define SR Ca concentration, in mM
    Casr = Casr_vec(i_sim);
    
    % Initial conditions for bond graph model
    X0 = [150e-6*W_i; ... Cai (fmol)
        Casr*W_sr; ... Casr (fmol)
        1e-4*W_isr; ... H (fmol)
        0.0363*W_i; ... MgADP (fmol)
        0.1*W_i; ... MgATP (fmol)
        1/9; ... P1 (fmol)
        15*W_i; ... Pi (fmol)
        1/9; ... P2 (fmol)
        1/9; ... P2a (fmol)
        1/9; ... P4 (fmol)
        1/9; ... P5 (fmol)
        1/9; ... P10 (fmol)
        1/9; ... P6 (fmol)
        1/9; ... P8 (fmol)
        1/9]; %P9 (fmol)
    
    % Simulate bond graph model to steady state, and record cycling rate
    [VOI, STATES, ALGEBRAIC] = SERCA_ode(tspan,X0);
    v_ss_vec(i_sim) = ALGEBRAIC(end,42);
    
    % Calculate cycling rate of Tran model
    v_ss_Tran(i_sim) = Tran_SERCA_model(Casr);
    
    % Calculate forward and reverse affinities of bond graph model
    Af = 2*ALGEBRAIC(end,1)+ALGEBRAIC(end,9);
    Ar = 2*ALGEBRAIC(end,3)+ALGEBRAIC(end,7)+ALGEBRAIC(end,13)+ALGEBRAIC(end,5);
    
    % Calculate free energy of the transporter at steady state
    DG_vec(i_sim) = Ar-Af;
    
    % Calculate affinities associated with transport and ATP hydrolysis
    A_Ca = 2*ALGEBRAIC(end,1) - 2*ALGEBRAIC(end,3);
    A_ATP = ALGEBRAIC(end,9) - (ALGEBRAIC(end,7)+ALGEBRAIC(end,13)+ALGEBRAIC(end,5));
    
    % Calculate the efficiency of the transporter
    if v_ss_vec(i_sim) >= 0
        eff_vec(i_sim) = -A_Ca/A_ATP;
    else
        eff_vec(i_sim) = -A_ATP/A_Ca;
    end
end

% Calculate equilibrium point
G0_MgATP = 11900; % Free energy of MgATP hydrolysis in J/mol
Casr_eq = 1e3*sqrt(exp(-G0_MgATP/R/T) * (150e-6)^2*0.1/(0.0363*15*1e-4));

% Plot steady-state cycling rate against SR calcium, comparing the bond
% graph and kinetic models
h1 = figure;
plot(Casr_vec,v_ss_Tran,'c','LineWidth',2);
hold on;
plot(Casr_vec,v_ss_vec,'k--','LineWidth',2);
xlim([0 5]);
ylim([-0.2 1]);
xlabel('[Ca^{2+}]_{sr} (mM)');
ylabel('Cycling rate (s^{-1})');
legend('Kinetic','Bond graph');
set(gca,'FontSize',20);

% Plot free energy of cycle for the bond graph model
h2 = figure;
plot(Casr_vec,1e-3*DG_vec,'k','LineWidth',2);
xlim([0 5]);
xlabel('[Ca^{2+}]_{sr} (mM)');
ylabel('\DeltaG (kJ/mol)');
set(gca,'FontSize',20);

% Plot steady-state power consumption for the bond graph model
h3 = figure;
plot(Casr_vec,-1e-3*DG_vec.*v_ss_vec,'k','LineWidth',2);
xlim([0 5]);
xlabel('[Ca^{2+}]_{sr} (mM)');
ylabel('Power (kJ/mol/s)');
set(gca,'FontSize',20);

% Plot steady-state efficiency for the bond graph model
h4 = figure;
plot(Casr_vec,100*eff_vec,'k','LineWidth',2);
xlim([0 5]);
xlabel('[Ca^{2+}]_{sr} (mM)');
ylabel('Efficiency (%)');
set(gca,'FontSize',20);

%% Save figures
if save_figures
    print_figure(h1,output_dir,'SERCA_v_ss');
    print_figure(h2,output_dir,'SERCA_DG');
    print_figure(h3,output_dir,'SERCA_power');
    print_figure(h4,output_dir,'SERCA_efficiency');
end