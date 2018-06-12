% This scripts compares the output of the updated model to data (Figures 7-8).
% The run_optimisation option may be set to true to run the fitting process
% again.

clear;
clc;
close all;

%% Script options
run_optimisation = false;
marker_size = 30;

%% Set up directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep 'NaK' filesep];
code_dir = [main_dir 'code' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];
dataset_func_dir = [code_dir 'dataset_functions' filesep];

addpath(dataset_func_dir);

%% Load Nakao and Gadsby voltage data
file_path = [data_dir 'nakao_gadsby_points.csv'];
struct_Nakao_Gadsby_V_data = load_Nakao_Gadsby_V(file_path);

%% Load Nakao and Gadsby voltage data with whole-cell currents
file_path = [data_dir 'Nakao_Gadsby_NaeV_currents.csv'];
struct_Nakao_Gadsby_I_data = load_Nakao_Gadsby_V_currents(file_path);

%% Load Nakao and Gadsby extracellular potassium data
file_path = [data_dir 'nakao_gadsby_Ko_activation_points.csv'];
struct_Nakao_Gadsby_Ke_data = load_Nakao_Gadsby_Ke(file_path);

%% Load Hansen intracellular Na data
file_path = [data_dir 'Hansen_Nai.csv'];
struct_Hansen_Nai_data = load_Hansen_Nai(file_path);

%% Load Friedrich ATP data
file_path = [data_dir 'Friedrich_ATP.csv'];
struct_Friedrich_MgATP_data = load_Friedrich_ATP(file_path);

%% Compile all data sources into a struct
struct_NaK_data = struct('Nakao_Gadsby_V_data',struct_Nakao_Gadsby_V_data,...
    'Nakao_Gadsby_I_data',struct_Nakao_Gadsby_I_data,...
    'Nakao_Gadsby_Ke_data',struct_Nakao_Gadsby_Ke_data,...
    'Hansen_Nai_data',struct_Hansen_Nai_data,...
    'Friedrich_MgATP_data',struct_Friedrich_MgATP_data);

%% Formulate minimisation problem for Nakao and Gadsby data
error_func = @(params_vec)NaK_total_error(params_vec,struct_NaK_data);

A = [];
b = [];
Aeq = [];
beq = [];
lb = [zeros(13,1); -Inf; 0];
ub = inf(15,1);

options_unc = optimoptions('fminunc','MaxFunEvals',10000);
options_ps = optimoptions('particleswarm','UseParallel',true);

%% Run parameter optimisation
if run_optimisation
    fmin = 1e10;
    tic
    for i_runs = 1:10
        i_runs
        [params_vec,fval,exitflag,output] = particleswarm(error_func,15,lb,ub,options_ps);
        [params_vec,fval,exitflag,output,grad,hessian] = fminunc(error_func,params_vec,options_unc);
        toc
        if fval < fmin
            params_vec_min = params_vec;
            fmin = fval
            output_min = output;
            hessian_min = hessian;
        end
    end
    params_vec = params_vec_min;
    fval = fmin;
    output = output_min;
    hessian = hessian_min;
    % Save parameters
    if strcmp(dataset,'predicted')
        save([storage_dir 'fitting_results_predicted.mat'],'params_vec','options_ps','options_unc','hessian');
    else
        save([storage_dir 'kinetic_fitting_results.mat'],'params_vec','options_ps','options_unc','hessian');
    end
else
    load([storage_dir 'kinetic_fitting_results.mat']);
end

%% Plot comparison to Nakao and Gadsby voltage data
% Constants
struct_input.T = 310; % Unit K

% Species concentrations (in mM)
struct_input.Ki = 0;
struct_input.Ke = 5.4;
struct_input.Nai = 50;
struct_input.Nae = 150;
struct_input.MgATP = 10;
struct_input.MgADP = 0;
struct_input.Pi_total = 0;
struct_input.H_conc = 10^(-4.4);

% Membrane potential (volts)
struct_input.V = -0.08;

% Compare simulated cycling rates to data (with normalisation)
h1 = figure;
hold on;
box on;
plot(1000*struct_Nakao_Gadsby_V_data.V_m,struct_Nakao_Gadsby_V_data.cyc_Nae1p5,'.c','MarkerSize',marker_size);
plot(1000*struct_Nakao_Gadsby_V_data.V_m,struct_Nakao_Gadsby_V_data.cyc_Nae50,'.r','MarkerSize',marker_size);
plot(1000*struct_Nakao_Gadsby_V_data.V_m,struct_Nakao_Gadsby_V_data.cyc_Nae100,'.b','MarkerSize',marker_size)
plot(1000*struct_Nakao_Gadsby_V_data.V_m,struct_Nakao_Gadsby_V_data.cyc_Nae150,'.g','MarkerSize',marker_size)

Nae_vec = [1.5; 50; 100; 150];
Vm_vec = (-120:1:60)/1000;
array_colour_str = {'c';'r';'b';'g'};
h_plot = zeros(1,length(Nae_vec));

array_normalising_factor = zeros(length(Nae_vec),1);
for i_Nae = 1:length(Nae_vec)
    Nae = Nae_vec(i_Nae);
    struct_input.Nae = Nae;
    struct_input.V = 0.04;
    array_normalising_factor(i_Nae) = 55/NaK_fitting_vcyc(params_vec,struct_input);
end

for i_Nae = 1:length(Nae_vec)
    Nae = Nae_vec(i_Nae);
    cycling_velocity_vec = zeros(length(Vm_vec),1);
    for i_Vm = 1:length(Vm_vec)
        V_m = Vm_vec(i_Vm);
        struct_input.Nae = Nae;
        struct_input.V = V_m;
        
        cycling_velocity_vec(i_Vm) = NaK_fitting_vcyc(params_vec,struct_input);
    end
    h_plot(1,i_Nae) = plot(1000*Vm_vec,cycling_velocity_vec*array_normalising_factor(i_Nae),[array_colour_str{i_Nae} '-'],'LineWidth',2);
end
xlabel('Membrane voltage (mV)');
ylabel('Cycling velocity (s^{-1})');
set(gca,'FontSize',18);
lgd = legend(h_plot,'[Na^+]_e = 1.5mM','[Na^+]_e = 50mM','[Na^+]_e = 100mM','[Na^+]_e = 150mM',...
    'Location','southeast');
lgd.FontSize = 12;
print_figure(h1,output_dir,'Terkildsen_kinetic_fit_KA_comparison');

% Compare simulated cycling rates to data (without normalisation)
h1b = figure;
hold on;
box on;
for i_Nae = 1:length(Nae_vec)
    Nae = Nae_vec(i_Nae);
    cycling_velocity_vec = zeros(length(Vm_vec),1);
    for i_Vm = 1:length(Vm_vec)
        V_m = Vm_vec(i_Vm);
        struct_input.Nae = Nae;
        struct_input.V = V_m;
        
        cycling_velocity_vec(i_Vm) = NaK_fitting_vcyc(params_vec,struct_input);
    end
    h_plot(1,i_Nae) = plot(1000*Vm_vec,cycling_velocity_vec,[array_colour_str{i_Nae} '-'],'LineWidth',2);
end
xlabel('Membrane voltage (mV)');
ylabel('Cycling velocity (s^{-1})');
set(gca,'FontSize',18);
lgd = legend(h_plot,'[Na^+]_e = 1.5mM','[Na^+]_e = 50mM','[Na^+]_e = 100mM','[Na^+]_e = 150mM',...
    'Location','southeast');
lgd.FontSize = 12;
print_figure(h1b,output_dir,'Terkildsen_vcyc');

%% Plot comparison to Nakao and Gadsby voltage data with whole-cell currents
struct_input.T = 310;  % Unit K

% Species concentrations (in mM)
struct_input.Ki = 0;
struct_input.Ke = 5.4;
struct_input.Nai = 50;
struct_input.Nae = 150;
struct_input.MgATP = 10;
struct_input.MgADP = 0;
struct_input.Pi_total = 0;
struct_input.H_conc = 10^(-4.4);

% Membrane potential (volts)
struct_input.V = -0.08;

Nae_vec = [50; 100; 150];

pump_density = params_vec(15);

array_colour_str = {'r';'b';'g'};
h_plot = zeros(1,length(Nae_vec));

h = figure;
hold on;
box on;
Vm_vec = struct_Nakao_Gadsby_I_data.V_m;
for i_Nae = 1:length(Nae_vec)
    plot(1000*Vm_vec,struct_Nakao_Gadsby_I_data.I(:,i_Nae),[array_colour_str{i_Nae} '.'],...
        'MarkerSize',marker_size);
end

Vm_vec = transpose((-140:1:60)/1000);
for i_Nae = 1:length(Nae_vec)
    Nae = Nae_vec(i_Nae);
    I_model_vec = zeros(length(Vm_vec),1);
    % Change Nae concentration
    struct_input.Nae = Nae;
    struct_input.V = 0.04;
    for i_Vm = 1:length(Vm_vec)
        V_m = Vm_vec(i_Vm);
        % Change membrane voltage
        struct_input.V = V_m;
        
        I_model = NaK_fitting_vcyc(params_vec,struct_input) * pump_density * 1.6e-5;
        I_model_vec(i_Vm) = I_model;
    end
    h_plot(1,i_Nae) = plot(1000*Vm_vec,I_model_vec,[array_colour_str{i_Nae} '-'],...
        'LineWidth',2);
end
xlabel('Membrane voltage (mV)');
ylabel('Current Density ({\mu}A/{\mu}F)');
% title('Comparison of fitted model to data');
set(gca,'FontSize',18);
lgd = legend(h_plot,'[Na^+]_e = 50mM','[Na^+]_e = 100mM','[Na^+]_e = 150mM',...
    'Location','southeast');
lgd.FontSize = 12;
print_figure(h,output_dir,'Terkildsen_fit_NG_I');

%% Plot comparison to Nakao and Gadsby extracellular potassium data
struct_input.T = 310;  % Unit K

% Species concentrations (in mM)
struct_input.Nae = 150;
struct_input.Nai = 50;
struct_input.Ki = 140;
struct_input.MgATP = 10;
struct_input.MgADP = 0.02;
struct_input.Pi_total = 0.5;
struct_input.H_conc = 10^(-4.4);

% Membrane potential (volts)
struct_input.V = 0;

% Model parameters 
Ke_vec = transpose(0.1:0.1:10);
num_data_points = length(Ke_vec);

struct_input.Ke = 5.4; % Unit mM
normalising_factor = NaK_fitting_vcyc(params_vec,struct_input);
cycling_velocity_model_vec = zeros(num_data_points,1);

for i_Ke = 1:num_data_points
    % Change Ke
    struct_input.Ke = Ke_vec(i_Ke);
    cycling_velocity_model_vec(i_Ke) = NaK_fitting_vcyc(params_vec,struct_input)/normalising_factor;
end

h = figure;
hold on;
box on;
plot(struct_Nakao_Gadsby_Ke_data.Ke,struct_Nakao_Gadsby_Ke_data.cyc_rate,...
    '.k','MarkerSize',marker_size);
plot(Ke_vec,cycling_velocity_model_vec,'k-','LineWidth',1.5);
xlabel('[K^+]_e (mM)');
ylabel('Relative cycling velocity');
% title('Comparison of fitted model to data');
set(gca,'FontSize',18);
print_figure(h,output_dir,'Terkildsen_fit_NG_Ke_comparison');

%% Plot comparison to Hansen intracellular sodium data
% Species amounts
% Constants
struct_input.T = 308; % Unit K

% Species concentrations (in mM)
struct_input.Nae = 0;
struct_input.Ke = 15;
struct_input.Ki = 80;
struct_input.MgATP = 2;
struct_input.MgADP = 0;
struct_input.Pi_total = 1;
struct_input.H_conc = 10^(-4.2);

% Membrane potential (volts)
struct_input.V = 0;

Nai_vec = transpose(1:1:80);
num_data_points = length(Nai_vec);

struct_input.Nai = 50;
normalising_factor_model = NaK_fitting_vcyc(params_vec,struct_input);
normalising_factor_data = struct_Hansen_Nai_data.cyc_rate(end-1);
cycling_velocity_model_vec = zeros(num_data_points,1);

for i_Nai = 1:num_data_points
    struct_input.Nai = Nai_vec(i_Nai);
    cycling_velocity_model_vec(i_Nai) = NaK_fitting_vcyc(params_vec,struct_input)/normalising_factor_model;
end

h = figure;
hold on;
box on;
plot(struct_Hansen_Nai_data.Nai,struct_Hansen_Nai_data.cyc_rate/normalising_factor_data,...
    '.k','MarkerSize',marker_size);
plot(Nai_vec,cycling_velocity_model_vec,'k-','LineWidth',1.5);
xlabel('[Na^+]_i (mM)');
ylabel('Relative cycling velocity');
% title('Comparison of fitted model to data');
set(gca,'FontSize',18);
print_figure(h,output_dir,'Terkildsen_fit_Hansen_Nai_comparison');

%% Plot comparison to Friedrich ATP data
% Constants
struct_input.T = 297; % Unit K

% Species concentrations (in mM)
struct_input.Nae = 0;
struct_input.Nai = 40;
struct_input.Ke = 5;
struct_input.Ki = 0;
struct_input.MgADP = 0;
struct_input.Pi_total = 0;
struct_input.H_conc = 10^(-4.4);

% Membrane potential (volts)
struct_input.V = 0;

MgATP_vec = 0.001:0.001:0.6;
num_data_points = length(MgATP_vec);

struct_input.MgATP = MgATP_vec(end);
normalising_factor = NaK_fitting_vcyc(params_vec,struct_input);

cycling_velocity_model_vec = zeros(num_data_points,1);

for i_MgATP = 1:num_data_points
    struct_input.MgATP = MgATP_vec(i_MgATP);
    cycling_velocity_model_vec(i_MgATP) = NaK_fitting_vcyc(params_vec,struct_input)/normalising_factor;
end

h = figure;
hold on;
box on;
plot(struct_Friedrich_MgATP_data.MgATP,struct_Friedrich_MgATP_data.cyc_rate,...
    '.k','MarkerSize',marker_size);
plot(MgATP_vec,cycling_velocity_model_vec,'k-','LineWidth',1.5);
xlabel('[MgATP] (mM)');
ylabel('Relative cycling velocity');
% title('Comparison of fitted model to data');
set(gca,'FontSize',18);
print_figure(h,output_dir,'Terkildsen_fit_Friedrich_MgATP_comparison');
