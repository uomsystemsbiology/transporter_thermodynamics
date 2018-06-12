% This script simulates the electrogenic transporter model described in the
% manuscript, and generates the plots in Figure 3.

clear;
clc;
close all;

%% Set directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];

%% Options
save_figures = true;

%% Define constants
R = 8.314; % Unit J/mol/K
T = 310; % Unit K
F = 96485; % Unit C/mol

%% Run system to steady state, for different voltages
V_vec = transpose((-100:1:100)/1000);
num_sim = length(V_vec);
v_ss_array = zeros(num_sim,1);
DG_vec = zeros(num_sim,1);

tspan = [0 1];

for i_sim = 1:num_sim
    % Set initial voltage and charge
    V = V_vec(i_sim);
    q0 = F^2*V;

    % Initial conditions
    X0 = [1; ...    % E1
        1; ...      % E2
        100; ...    % Se
        10; ...     % Si
        q0];        % q_mem

    % Run simulation to steady state, and record cycling rate and free
    % energy
    [VOI, STATES, ALGEBRAIC] = charged_species_ode(tspan,X0);
    v_ss_array(i_sim) = ALGEBRAIC(end,23)/(sum(X0(1:2)));
    DG_vec(i_sim) = ALGEBRAIC(1,5)-ALGEBRAIC(1,7)-F*ALGEBRAIC(1,9);
end

% Calculate the predicted equilibrium potential
V_eq = R*T/F * log(100/10);

% Plot cycling rate against membrane potential
h1 = figure;
plot(V_vec,v_ss_array,'k','LineWidth',2);
hold on;
plot([-0.1 0.1],[0 0],'b--','LineWidth',2);
plot(V_eq,0,'b.','MarkerSize',30);
xlabel('V_m');
ylabel('v_{cyc}');
set(gca,'FontSize',20);

% Plot free energy against membrane potential
h2 = figure;
plot(V_vec,DG_vec,'k','LineWidth',2);
hold on;
plot([-0.1 0.1],[0 0],'b--','LineWidth',2);
plot(V_eq,0,'b.','MarkerSize',30);
xlabel('V_m');
ylabel('\DeltaG');
set(gca,'FontSize',20);

% Plot free energy against cycling rate
h3 = figure;
plot(DG_vec,v_ss_array,'k','LineWidth',2);
hold on;
plot(0,0,'k.','MarkerSize',30);
text(1000,0,'\DeltaG = 0, v_{ss} = 0','FontSize',16,...
    'HorizontalAlignment','left','VerticalAlignment','top');
xlabel('\DeltaG');
ylabel('v_{cyc}');
set(gca,'FontSize',20);


%% Save figures
if save_figures
    print_figure(h1,output_dir,'charged_species_V');
    print_figure(h2,output_dir,'charged_species_DG');
    print_figure(h3,output_dir,'charged_species_DG_vss');
end