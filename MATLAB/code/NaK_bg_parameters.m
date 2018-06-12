% This script loads the kinetic parameters from fitting the Na/K ATPase
% kinetic model, and solves for the bond graph parameters.

clear;
clc;
close all;

%% Define volumes using data from Luo and Rudy (unit pL)
W_i = 38; % Intracelluar volume
W_e = 5.182; % Extracellular volume

%% Directory names
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep 'NaK' filesep];
output_dir = [main_dir 'output' filesep];
storage_dir = [main_dir 'storage' filesep];

%% Load forward matrix
stoich_path = [data_dir 'forward_matrix.txt'];
stoich_file_id = fopen(stoich_path,'r');

stoich_file_data = textscan(stoich_file_id,'%s','delimiter','\n');
fclose(stoich_file_id);

num_rows = length(stoich_file_data{1});
num_cols = sum(stoich_file_data{1}{1} == ',')+1;

N_f = zeros(num_rows,num_cols);

for i_row = 1:num_rows
    line_str = stoich_file_data{1}{i_row};
    line_split = regexp(line_str,',','split');
    
    for i_col = 1:num_cols
        N_f(i_row,i_col) = str2double(line_split{i_col});
    end
end

N_fT = transpose(N_f);

%% Load reverse matrix
stoich_path = [data_dir 'reverse_matrix.txt'];
stoich_file_id = fopen(stoich_path,'r');

stoich_file_data = textscan(stoich_file_id,'%s','delimiter','\n');
fclose(stoich_file_id);

num_rows = length(stoich_file_data{1});
num_cols = sum(stoich_file_data{1}{1} == ',')+1;

N_r = zeros(num_rows,num_cols);

for i_row = 1:num_rows
    line_str = stoich_file_data{1}{i_row};
    line_split = regexp(line_str,',','split');
    
    for i_col = 1:num_cols
        N_r(i_row,i_col) = str2double(line_split{i_col});
    end
end

N_rT = transpose(N_r);

%% Calculate stoichiometric matrix
N = N_r - N_f;
N_T = N_rT - N_fT;

I = eye(num_cols);

M = [I N_fT; I N_rT];
M_rref = rref(M);
M_rank = rank(M);

M_T = transpose(M);
M_T_rref = rref(M_T);

%% Calculate the bond graph constants 
% Load kinetic parameters
load([storage_dir 'kinetic_fitting_results.mat']);
params = params_vec;

% Use dissociation constants to calculate fast kinetic constants
fast_kinetic_constant = 1e6;

K_d_Nai = params(3);
k_3b_p = fast_kinetic_constant;
k_3b_m = fast_kinetic_constant*(0.5*K_d_Nai);

k_4b_p = fast_kinetic_constant;
k_4b_m = fast_kinetic_constant*(2*K_d_Nai);

K_d_Nai_0 = params(1);
k_5b_m_0 = fast_kinetic_constant;
k_5b_p_0 = fast_kinetic_constant/K_d_Nai_0;

k_6b_p = params(7);
k_6b_m = k_6b_p/6.3;

k_7b_p = params(8);
k_7b_m = params(11);

K_d_Nae_0 = params(2);
k_8b_m_0 = fast_kinetic_constant;
k_8b_p_0 = fast_kinetic_constant*K_d_Nae_0;

K_d_Nae = params(4);
k_9b_m = fast_kinetic_constant;
k_9b_p = fast_kinetic_constant*2*K_d_Nae;

k_10b_m = fast_kinetic_constant;
k_10b_p = fast_kinetic_constant*0.5*K_d_Nae;

K_d_Ke = params(5);
k_11b_p = fast_kinetic_constant;
k_11b_m = fast_kinetic_constant*0.5*K_d_Ke;

k_12b_p = fast_kinetic_constant;
k_12b_m = fast_kinetic_constant*2*K_d_Ke;

k_13b_p = params(9);
k_13b_m = params(12);

K_d_MgATP = params(6);
k_14b_p = fast_kinetic_constant;
k_14b_m = fast_kinetic_constant*K_d_MgATP;

k_15b_p = params(10);
k_15b_m = params(13);

% Calculate remaining parameter
G_MgATP_0 = 11900;
RT = 8.314*310;
K_MgATP = exp(-G_MgATP_0/RT)*10^6;
K_d_Ki = sqrt((k_3b_m*k_4b_m*k_5b_m_0*k_6b_m*k_7b_m*k_8b_m_0*k_9b_m*k_10b_m*k_11b_m*k_12b_m*k_13b_m*k_14b_m*k_15b_m*K_MgATP)/(k_3b_p*k_4b_p*k_5b_p_0*k_6b_p*k_7b_p*k_8b_p_0*k_9b_p*k_10b_p*k_11b_p*k_12b_p*k_13b_p*k_14b_p*k_15b_p));
% K_d_Ki = 255.13;

k_1b_m = fast_kinetic_constant;
k_1b_p = fast_kinetic_constant*2*K_d_Ki;

k_2b_m = fast_kinetic_constant;
k_2b_p = fast_kinetic_constant*0.5*K_d_Ki;

% Define k
k = transpose([k_1b_p k_2b_p k_3b_p k_4b_p k_5b_p_0 k_6b_p k_7b_p k_8b_p_0 ...
    k_9b_p k_10b_p k_11b_p k_12b_p k_13b_p k_14b_p k_15b_p ...
    k_1b_m k_2b_m k_3b_m k_4b_m k_5b_m_0 k_6b_m k_7b_m k_8b_m_0 ...
    k_9b_m k_10b_m k_11b_m k_12b_m k_13b_m k_14b_m k_15b_m]);

% Check if the thermodynamic constraint holds
K_MgATP_est = prod(k(1:15))/prod(k(16:30));
G_MgATP_0_est = -RT*log(K_MgATP_est/1e6);

%% Add extra constraints for sodium, potassium and MgATP hydrolysis
K_Na = 1;
M_Na = zeros(1,size(M,2));
M_Na(33) = 1;
M_Na(34) = -1;

K_K = 1;
M_K = zeros(1,size(M,2));
M_K(31) = 1;
M_K(32) = -1;

M_MgATP = zeros(1,size(M,2));
M_MgATP(35) = 1;
M_MgATP(36) = -1;
M_MgATP(37) = -1;
M_MgATP(38) = -1;

k_mod = [k; K_Na; K_K; K_MgATP];
M_mod = [M; M_Na; M_K; M_MgATP];

%% Calculate bond graph constants
% Note: units of kappa are mmol/s, units of K are mmol^-1
lambdaW = exp(pinv(M_mod) * log(k_mod));
W = [ones(30,1); W_i; W_e; W_i; W_e; W_i; W_i; W_i; W_i];

lambda = lambdaW./W;

% Check that equation has been solved correctly
k_check = exp(M*log(lambdaW));
kappa = lambda(1:num_cols);
K = lambda(num_cols+1:end);

diff = sum(abs((k_mod - exp(M_mod*log(lambdaW)))./k_mod));

%% Save results
save([storage_dir 'NaK_bg_params.mat'],'kappa','K','params');
