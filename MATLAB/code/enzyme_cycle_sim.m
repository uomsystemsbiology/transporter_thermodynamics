% This script simulates the enzyme cycle model described in the
% manuscript, and generates the plots in Figure 3.

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

%% Plot the transients of enzyme cycle
% Time span
tspan = [0 0.2];

% Initial conditions
X0 = [1; ...    % E1
    1; ...      % E2
    10; ...     % Se
    100];       % Si

% Run simulation
[VOI, STATES, ALGEBRAIC, CONSTANTS] = enzyme_cycle_ode(tspan,X0);

t = VOI;
x_E1 = STATES(:,1);
x_E2 = STATES(:,2);

% Plot transporter states against time
h1 = figure;
plot(t,x_E1,t,x_E2,'LineWidth',2);
xlabel('Time (s)');
ylabel('Molar amount');
legend('x_{E1}','x_{E2}');
set(gca,'FontSize',20);

v1 = ALGEBRAIC(:,19);
v2 = ALGEBRAIC(:,28);

% Plot reaction velocities against time
h2 = figure;
plot(t,v1,t,v2,'LineWidth',2);
hold on;
plot([0 0.2],[v1(end) v1(end)],'k--','LineWidth',2);
ylim([-2 4]);
xlabel('Time (s)');
ylabel('Reaction velocity');
legend('v_1','v_2');
set(gca,'FontSize',20);

%% Plot steady-state cycling velocities against Se
Se_vec = transpose(1:200);
num_sim = length(Se_vec);
v_ss_vec = zeros(num_sim,1);

% Time span
tspan = [0 0.4];

for i_sim = 1:num_sim
    % Initial conditions
    Se = Se_vec(i_sim);
    X0 = [1; ...    % E1
        1; ...      % E2
        Se; ...     % Se
        100];       % Si
    
    % Run simulation, and record steady-state cycling rate
    [VOI, STATES, ALGEBRAIC] = enzyme_cycle_ode(tspan,X0);
    v_ss_vec(i_sim) = ALGEBRAIC(end,19)/(X0(1)+X0(2));
end

% Plot cycling rate against Se
h3 = figure;
plot(Se_vec,v_ss_vec,'k','LineWidth',2);
hold on;
plot([0 100],[0 0],'b--','LineWidth',2);
plot([100 100],[-0.5 0],'b--','LineWidth',2);
plot(100,0,'b.','MarkerSize',30);
text(100,0,'x_{Si} = x_{Se}','FontSize',20,...
    'HorizontalAlignment','left','VerticalAlignment','bottom');
xlabel('x_{Se}');
ylabel('v_{cyc}');
set(gca,'FontSize',20);

%% Save figures
if save_figures
    print_figure(h1,output_dir,'enzyme_cycle_x');
    print_figure(h2,output_dir,'enzyme_cycle_v');
    print_figure(h3,output_dir,'enzyme_cycle_ss');
end