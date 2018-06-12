% This script simulates the coupled transporter model described in the
% manuscript, and generates the plots in Figure 2.

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

%% Run system to steady state, for different values of A
Se_vec = transpose(1:200);
A_vec = [2; 5; 10; 20];
num_sim = length(Se_vec);
num_A = length(A_vec);
v_ss_array = zeros(num_sim,num_A);

tspan = [0 5];

for i_A = 1:num_A
    A = A_vec(i_A);
    for i_sim = 1:num_sim
        Se = Se_vec(i_sim);

        % Initial conditions
        X0 = [A; ...    % A
            1; ...      % B
            0.5; ...     % E1
            0.5; ...       % E2
            0.5; ...    % E3
            0.5; ...      % E4
            Se; ...     % Se
            10];       % Si
        
        % Run simulation, and record steady-state cycling rates
        [VOI, STATES, ALGEBRAIC] = coupled_reactions_ode(tspan,X0);
        v_ss_array(i_sim,i_A) = ALGEBRAIC(end,35)/(sum(X0(3:6)));
    end
end

% Plot cycling rate against Se for different values of A
h1 = figure;
h = plot(Se_vec,v_ss_array,'LineWidth',2);
h(1).Color = 'c';
h(2).Color = 'r';
h(3).Color = 'b';
h(4).Color = 'g';
hold on;
plot([0 200],[0 0],'k--','LineWidth',2);
plot(20,0,'c.','MarkerSize',30);
plot(50,0,'r.','MarkerSize',30);
plot(100,0,'b.','MarkerSize',30);
plot(200,0,'g.','MarkerSize',30);

ylim([-0.05 0.10]);
xlabel('x_{Se}');
ylabel('v_{cyc}');
legend('x_A = 2','x_A = 5','x_A = 10','x_A = 20');
set(gca,'FontSize',20);

%% Run system to steady state, examining effects of A
A_vec = transpose(1:100);
num_A = length(A_vec);
v_ss_array = zeros(num_A,1);

tspan = [0 5];

for i_A = 1:num_A
    A = A_vec(i_A);
        
    % Initial conditions
    X0 = [A; ...    % A
        1; ...      % B
        0.5; ...     % E1
        0.5; ...       % E2
        0.5; ...    % E3
        0.5; ...      % E4
        100; ...     % Se
        10];       % Si
    
    % Run simulation, and record steady-state cycling rates
    [VOI, STATES, ALGEBRAIC] = coupled_reactions_ode(tspan,X0);
    v_ss_array(i_A) = ALGEBRAIC(end,35)/(sum(X0(3:6)));
end


% Plot cycling rate against A
h2 = figure;
plot(A_vec,v_ss_array,'k','LineWidth',2);
hold on;
plot([0 10 10],[0 0 -0.04],'b--','LineWidth',2);
plot(10,0,'b.','MarkerSize',30);
text(13,0,'x_{A} = 10','FontSize',20,...
    'HorizontalAlignment','left','VerticalAlignment','top');
xlabel('x_A');
ylabel('v_{cyc}');
set(gca,'FontSize',20);

%% Examine relationship between free energy and steady state cycling rate

% Varying A
A_vec = transpose(1:100);
num_A = length(A_vec);
DG_array_A = zeros(num_A,1);
v_ss_array_A = zeros(num_A,1);

tspan = [0 5];

for i_A = 1:num_A
    A = A_vec(i_A);
        
    X0 = [A; ...    % A
        1; ...      % B
        0.5; ...     % E1
        0.5; ...       % E2
        0.5; ...    % E3
        0.5; ...      % E4
        100; ...     % Se
        10];       % Si
    [VOI, STATES, ALGEBRAIC] = coupled_reactions_ode(tspan,X0);

    DG_array_A(i_A) = ALGEBRAIC(1,13)+ALGEBRAIC(1,3)-ALGEBRAIC(1,15)-ALGEBRAIC(1,1);
    v_ss_array_A(i_A) = ALGEBRAIC(end,35)/(sum(X0(3:6)));
end

% Varying B
B_vec = transpose(1:100);
num_B = length(B_vec);
DG_array_B = zeros(num_B,1);
v_ss_array_B = zeros(num_B,1);

for i_B = 1:num_B
    B = A_vec(i_B);
        
    X0 = [100; ...    % A
        B; ...      % B
        0.5; ...     % E1
        0.5; ...       % E2
        0.5; ...    % E3
        0.5; ...      % E4
        100; ...     % Se
        10];       % Si
    [VOI, STATES, ALGEBRAIC] = coupled_reactions_ode(tspan,X0);

    DG_array_B(i_B) = ALGEBRAIC(1,13)+ALGEBRAIC(1,3)-ALGEBRAIC(1,15)-ALGEBRAIC(1,1);
    v_ss_array_B(i_B) = ALGEBRAIC(end,35)/(sum(X0(3:6)));
end

% Varying Se
Se_vec = transpose(logspace(2,4,101));
num_Se = length(Se_vec);
DG_array_Se = zeros(num_Se,1);
v_ss_array_Se = zeros(num_Se,1);

for i_Se = 1:num_Se
    Se = Se_vec(i_Se);
        
    X0 = [100; ...    % A
        1; ...      % B
        0.5; ...     % E1
        0.5; ...       % E2
        0.5; ...    % E3
        0.5; ...      % E4
        Se; ...     % Se
        10];       % Si
    [VOI, STATES, ALGEBRAIC] = coupled_reactions_ode(tspan,X0);

    DG_array_Se(i_Se) = ALGEBRAIC(1,13)+ALGEBRAIC(1,3)-ALGEBRAIC(1,15)-ALGEBRAIC(1,1);
    v_ss_array_Se(i_Se) = ALGEBRAIC(end,35)/(sum(X0(3:6)));
end

% Plot free energy against cycling rate, when changing A, B and Se
h3 = figure;
plot(DG_array_A,v_ss_array_A,'r','LineWidth',2);
hold on;
plot(DG_array_B,v_ss_array_B,'g','LineWidth',2);
plot(DG_array_Se,v_ss_array_Se,'b','LineWidth',2);
plot(0,0,'k.','MarkerSize',30);
text(0,0,'\DeltaG = 0, v_{ss} = 0','FontSize',16,...
    'HorizontalAlignment','right','VerticalAlignment','top');
ylim([-0.01 0.01]);
xlabel('\DeltaG');
ylabel('v_{cyc}');
legend('Varying A','Varying B','Varying Se');
set(gca,'FontSize',20);

%% Save figures
if save_figures
    print_figure(h1,output_dir,'coupled_reactions_Se');
    print_figure(h2,output_dir,'coupled_reactions_A');
    print_figure(h3,output_dir,'coupled_reactions_DG');
end