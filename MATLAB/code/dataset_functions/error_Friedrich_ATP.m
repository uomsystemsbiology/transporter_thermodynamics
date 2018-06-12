function final_error = error_Friedrich_ATP(params_vec,struct_Friedrich_MgATP_data)
% Function that calculates the error of a parameter set to the
% MgATP data of Friedrich et al.

% Constants
struct_input.T = 297;

% Species concentrations
struct_input.Nae = 0;
struct_input.Nai = 40;
struct_input.Ke = 5;
struct_input.Ki = 0;
struct_input.MgADP = 0;
struct_input.Pi_total = 0;
struct_input.H_conc = 10^(-4.4);
struct_input.V = 0;

MgATP_vec = struct_Friedrich_MgATP_data.MgATP;
num_data_points = length(MgATP_vec);

struct_input.MgATP = MgATP_vec(end);
normalising_factor = NaK_fitting_vcyc(params_vec,struct_input);

total_error = 0;

for i_MgATP = 1:num_data_points
    struct_input.MgATP = MgATP_vec(i_MgATP);

    weight = 1;

    cycling_velocity_data = struct_Friedrich_MgATP_data.cyc_rate(i_MgATP);
    cycling_velocity_model = NaK_fitting_vcyc(params_vec,struct_input)/normalising_factor;

    error = (cycling_velocity_data-cycling_velocity_model)^2;
    total_error = total_error + weight*error;
end
final_error = total_error;

end