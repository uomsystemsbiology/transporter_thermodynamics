% This script simulates the bond graph version of the updated Terkildsen et
% al. (2007) model, and generates the plots in Figure 9 of the manuscript.

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
T = 310; % Unit K
F = 96485; % Unit C/mol

W_i = 38; % Intracellular volume, unit pL
W_e = 5.182; % Extracellular volume, unit pL

A_cap = 1.534e-4; % Membrane area, unit cm^2
mem_cap = 1e9 * A_cap; % Membrane capacitance, unit fF

%% Run system to steady state for different Nae, V
V_vec = transpose((-120:10:60)/1000);
Nae_vec = transpose([1.5; 50; 100; 150]);
num_V = length(V_vec);
num_Nae = length(Nae_vec);

v_ss_vec = zeros(num_V,num_Nae);
v_ss_kinetic = zeros(num_V,num_Nae);

tspan = [0 1];

for i_Nae = 1:num_Nae
    Nae = Nae_vec(i_Nae);
    for i_V = 1:num_V
        V = V_vec(i_V); % Membrane potential, unit volts
                
        % Initial conditions
        X0 = [10^(-4.4)*W_i; ... H (fmol)
            5.4*W_e; ... Ke (fmol)
            0.001*W_i; ... Ki (fmol)
            0.001*W_i; ... MgADP (fmol)
            10*W_i; ... MgATP (fmol)
            Nae*W_e; ... Nae (fmol)
            50*W_i; ... Nai (fmol)
            1/15; ... P7 (fmol)
            0.001*W_i; ... Pi (fmol)
            V*mem_cap; ... mem (fC)
            1/15; ... P14 (fmol)
            1/15; ... P15 (fmol)
            1/15; ... P1 (fmol)
            1/15; ... P2 (fmol)
            1/15; ... P3 (fmol)
            1/15; ... P4 (fmol)
            1/15; ... P5 (fmol)
            1/15; ... P6 (fmol)
            1/15; ... P10 (fmol)
            1/15; ... P11 (fmol)
            1/15; ... P12 (fmol)
            1/15; ... P13 (fmol)
            1/15; ... P8 (fmol)
            1/15]; %P9 (fmol)
        
        % Simulate bond graph model to steady state, and record cycing rate
        [VOI, STATES, ALGEBRAIC] = NaK_ode(tspan,X0);
        v_ss_vec(i_V,i_Nae) = ALGEBRAIC(end,130);
        
        % Calculate steady-state cycling rate for kinetic model
        v_ss_kinetic(i_V,i_Nae) = NaK_kinetic_ode(1000*V,Nae);
    end
end
    
% Compare the steady-state cycling rates of the kinetic and bond graph
% models
h1 = figure;
h_kinetic = plot(1000*V_vec,v_ss_kinetic,'c','LineWidth',2);
hold on;
h_bg = plot(1000*V_vec,v_ss_vec,'k--','LineWidth',2);
xlim([-160 90]);
xlabel('V_m (mV)');
ylabel('Cycling rate (s^{-1})');
legend([h_kinetic(1) h_bg(1)],'Kinetic','Bond graph','Location','southeast');
set(gca,'FontSize',20);

%% Check the equilibrium voltages of the pump
V_vec = transpose((-400:10:-200)/1000);
MgATP_vec = transpose([1 6.95]);
num_V = length(V_vec);
num_MgATP = length(MgATP_vec);

G_MgATP = 11900; % Free energy of MgATP hydrolysis in J/mol

v_ss_vec = zeros(num_V,num_MgATP);
V_eq = zeros(1,num_MgATP);

tspan = [0 1];

for i_MgATP = 1:num_MgATP
    MgATP = MgATP_vec(i_MgATP);
    
    % Calculate the predicted equilibrium potential
    V_eq(i_MgATP) = 1/F * G_MgATP ...
        + R*T/F *(log(0.035*0.3971*10^(-4.095)/MgATP*1e-6) + 3*log(140/10) + 2*log(145/5.4));
    for i_V = 1:num_V
        V = V_vec(i_V); % Membrane potential, unit volts
        
        % Initial conditions
        X0 = [10^(-4.095)*W_i; ... H (fmol)
            5.4*W_e; ... Ke (fmol)
            145*W_i; ... Ki (fmol)
            0.035*W_i; ... MgADP (fmol)
            MgATP*W_i; ... MgATP (fmol)
            140*W_e; ... Nae (fmol)
            10*W_i; ... Nai (fmol)
            1/15; ... P7 (fmol)
            0.3971*W_i; ... Pi (fmol)
            V*mem_cap; ... mem (fC)
            1/15; ... P14 (fmol)
            1/15; ... P15 (fmol)
            1/15; ... P1 (fmol)
            1/15; ... P2 (fmol)
            1/15; ... P3 (fmol)
            1/15; ... P4 (fmol)
            1/15; ... P5 (fmol)
            1/15; ... P6 (fmol)
            1/15; ... P10 (fmol)
            1/15; ... P11 (fmol)
            1/15; ... P12 (fmol)
            1/15; ... P13 (fmol)
            1/15; ... P8 (fmol)
            1/15]; %P9 (fmol)
        
        % Simulate bond graph model to steady state, and record cycing rate
        [VOI, STATES, ALGEBRAIC] = NaK_ode(tspan,X0);
        v_ss_vec(i_V,i_MgATP) = ALGEBRAIC(end,130);

    end
end
    
% Plot steady-state cycling rate of bond graph model, with predicted
% equilibrium points
h2 = figure;
h = plot(1000*V_vec,v_ss_vec,'LineWidth',2);
hold on;
h(1).Color = 'b';
h(2).Color = 'r';
plot([-400 -200],[0 0],'k--','LineWidth',2);
plot(1000*V_eq(1),0,'b.','MarkerSize',30);
plot(1000*V_eq(2),0,'r.','MarkerSize',30);
xlabel('V_m (mV)');
ylabel('Cycling rate (s^{-1})');
legend('[MgATP]=1mM','[MgATP]=6.95mM','Location','northwest');
set(gca,'FontSize',20);

%% Plot standard free energy against cycling rate, while varying MgATP and membrane potential
V_vec = transpose((-400:10:-200)/1000);
MgATP_vec = transpose(1:5);
num_V = length(V_vec);
num_MgATP = length(MgATP_vec);

v_ss_vec = zeros(num_V,num_MgATP);
DG_array = zeros(num_V,num_MgATP);

tspan = [0 1];

for i_MgATP = 1:num_MgATP
    MgATP = MgATP_vec(i_MgATP);
    for i_V = 1:num_V
        V = V_vec(i_V); % Membrane potential, unit volts
        
        % Initial conditions
        X0 = [10^(-4.095)*W_i; ... H (fmol)
            5.4*W_e; ... Ke (fmol)
            145*W_i; ... Ki (fmol)
            0.035*W_i; ... MgADP (fmol)
            MgATP*W_i; ... MgATP (fmol)
            140*W_e; ... Nae (fmol)
            10*W_i; ... Nai (fmol)
            1/15; ... P7 (fmol)
            0.3971*W_i; ... Pi (fmol)
            V*mem_cap; ... mem (fC)
            1/15; ... P14 (fmol)
            1/15; ... P15 (fmol)
            1/15; ... P1 (fmol)
            1/15; ... P2 (fmol)
            1/15; ... P3 (fmol)
            1/15; ... P4 (fmol)
            1/15; ... P5 (fmol)
            1/15; ... P6 (fmol)
            1/15; ... P10 (fmol)
            1/15; ... P11 (fmol)
            1/15; ... P12 (fmol)
            1/15; ... P13 (fmol)
            1/15; ... P8 (fmol)
            1/15]; %P9 (fmol)
        
        % Simulate bond graph model to steady state, and record cycing rate
        [VOI, STATES, ALGEBRAIC] = NaK_ode(tspan,X0);
        v_ss_vec(i_V,i_MgATP) = ALGEBRAIC(end,130);
        
        % Calculate the forward and reverse affinities of the chemical
        % components of the pump
        Af = ALGEBRAIC(end,17) + 3*ALGEBRAIC(end,21) + 2*ALGEBRAIC(end,11);
        Ar = ALGEBRAIC(end,15) + ALGEBRAIC(end,25) + ALGEBRAIC(end,9) + 3*ALGEBRAIC(end,19) + 2*ALGEBRAIC(end,13);
        
        % Calculate the free energy of the pump
        DG_array(i_V,i_MgATP) = Ar-Af - F*V;
    end
end

% Plot steady-state cycling rate against free energy
h3 = figure;
h = plot(1e-3*DG_array,v_ss_vec,'k','LineWidth',2);
xlabel('\DeltaG (kJ/mol)');
ylabel('Cycling rate (s^{-1})');
set(gca,'FontSize',20);

%% Save figures
if save_figures
    print_figure(h1,output_dir,'NaK_v_ss');
    print_figure(h2,output_dir,'NaK_eq');
    print_figure(h3,output_dir,'NaK_DG');
end