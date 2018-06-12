function final_error = error_Hansen_Nai(params_vec,struct_Hansen_Nai_data)
% Function that calculates the error of a parameter set to the
% intracellular Na data of Hansen

% Constants
struct_input.T = 308;


% Species concentrations
struct_input.Nae = 0;
struct_input.Ke = 15;
struct_input.Ki = 80;
struct_input.MgATP = 2;
struct_input.MgADP = 0;
struct_input.Pi_total = 1;
struct_input.H_conc = 10^(-4.2);
struct_input.V = 0;

Nai_vec = struct_Hansen_Nai_data.Nai;
num_data_points = length(Nai_vec);

struct_input.Nai = 50;
normalising_factor_model = NaK_fitting_vcyc(params_vec,struct_input);
normalising_factor_data = struct_Hansen_Nai_data.cyc_rate(end-1);

total_error = 0;

for i_Nai = 1:num_data_points
    struct_input.Nai = Nai_vec(i_Nai);

    weight = 1;
    if Nai_vec(i_Nai) <= 15
        weight = weight*2;
    end

    cycling_velocity_data = struct_Hansen_Nai_data.cyc_rate(i_Nai)/normalising_factor_data;
    cycling_velocity_model = NaK_fitting_vcyc(params_vec,struct_input)/normalising_factor_model;

    error = (cycling_velocity_data-cycling_velocity_model)^2;
    total_error = total_error + weight*error;
end
final_error = total_error;

end