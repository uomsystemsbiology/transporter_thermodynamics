% This script calculates the bond graph parameters for the SERCA model of
% Tran et al. (2009) by using the kinetic parameters and stoichiometric
% matrix.

clear;
clc;
close all;

%% Set directories
current_dir = cd;
Idx_backslash = find(current_dir == filesep);
main_dir = current_dir(1:Idx_backslash(end));
data_dir = [main_dir 'data' filesep 'SERCA' filesep];
output_dir = [main_dir 'output' filesep];

%% Define the volumes of each compartment (units pL)
W_i = 38.0; % Intracellular volume
W_sr = 2.28; % SR volume
W_isr = W_i + W_sr; % Intracellular volume + SR volume

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

N_cT = zeros(2,size(M,2));

% Add constraint for MgATP hydrolysis
N_cT(1,22) = 1;
N_cT(1,23) = -1;
N_cT(1,24) = -1;
N_cT(1,19) = -1;

% Add constraint for calcium transport
N_cT(2,20) = 1;
N_cT(2,21) = -1;

M = [M; N_cT];

%% Set the kinetic rate constants

% R1-2
k_12p = 25900;

% R2-4
K_dCai = 0.91;
k_24m = 259e8;
k_24p = k_24m/K_dCai^2;

% R2-2a
K_dH1 = 1.09e-5;
k_22am = 259e8;
k_22ap = k_22am/K_dH1;

% R4-5
K_dHi = 3.54e-3;
k_45p = 259e8;
k_45m = k_45p/K_dHi;

% R5-6
k_56p = 2540;
k_56m = 67200;

% R6-8
K_dCasr = 2.24;
k_68m = 259e8;
k_68p = k_68m*K_dCasr^2;

% R8-9
K_dHsr = 1.05e-8;
k_89m = 259e8;
k_89p = k_89m/K_dHsr;

% R9-10
K_dH = 7.24e-5;
k_910p = 259e8;
k_910m = k_910p/K_dH;

% R10-1
k_101p = 20.5;
k_101m = 149;

% 
G_0_MgATP = 11900;
R = 8.314;
T = 310;
K_MgATP = exp(-G_0_MgATP/(R*T))*10^6;

% Calculate remaining parameter using detailed balance
k_12m = (k_12p*k_24p*k_45p*k_56p*k_68p*k_89p*k_910p*k_101p)/(k_24m*k_45m*k_56m*k_68m*k_89m*k_910m*k_101m*K_MgATP);

%% Calculate bond graph constants from kinetic parameters
% Note: units of kappa are fmol/s, units of K are fmol^-1
k = transpose([k_12p k_24p k_22ap k_45p k_56p k_68p k_89p k_910p k_101p k_12m k_24m k_22am k_45m k_56m k_68m k_89m k_910m k_101m K_MgATP 1]);
lambdaW = exp(pinv(M) * log(k));

% Check that kinetic parameters are reproduced by bond graph parameters
k_sub = exp(M*log(lambdaW));
diff = (k - k_sub)./k;
error = sum(abs(diff));

% Check that there is a detailed balance constraint
Z = transpose(null(transpose(M),'r'));

%% Save bond graph parameters
W = [ones(18,1); W_isr; W_i; W_sr; W_i; W_i; W_i];
lambda = lambdaW./W;

kappa = lambda(1:num_cols);
K = lambda(num_cols+1:end);

save([output_dir 'SERCA_params.mat'],'kappa','K');
