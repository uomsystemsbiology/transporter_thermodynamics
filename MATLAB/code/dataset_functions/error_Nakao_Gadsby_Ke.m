function final_error = error_Nakao_Gadsby_Ke(params_vec,struct_Nakao_Gadsby_Ke_data)
% Function that calculates the error of a parameter set to the
% extracellular data of Nakao and Gadsby

struct_input.T = 310;

% Species concentrations
struct_input.Nae = 150;
struct_input.Nai = 50;
struct_input.Ki = 140;
struct_input.MgATP = 10;
struct_input.MgADP = 0.02;
struct_input.Pi_total = 0.5;
struct_input.H_conc = 10^(-4.4);
struct_input.V = 0;

Ke_vec = struct_Nakao_Gadsby_Ke_data.Ke;
num_data_points = length(Ke_vec);

struct_input.Ke = 5.4;
normalising_factor = NaK_fitting_vcyc(params_vec,struct_input);

total_error = 0;

for i_Ke = 1:num_data_points
    struct_input.Ke = Ke_vec(i_Ke);

    weight = 1;
    if Ke_vec(i_Ke) > 4.5
        weight = weight*15;
    end

    cycling_velocity_data = struct_Nakao_Gadsby_Ke_data.cyc_rate(i_Ke);
    cycling_velocity_model = NaK_fitting_vcyc(params_vec,struct_input)/normalising_factor;

    error = (cycling_velocity_data-cycling_velocity_model)^2;
    total_error = total_error + weight*error;
end
final_error = total_error;

end