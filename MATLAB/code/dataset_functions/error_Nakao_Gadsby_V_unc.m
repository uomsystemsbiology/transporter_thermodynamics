function final_error = error_Nakao_Gadsby_V_unc(params_vec,struct_Nakao_Gadsby_V_data)

% Constants
struct_input.T = 310;

% Species concentrations
struct_input.Ki = 0;
struct_input.Ke = 5.4;
struct_input.Nai = 50;
struct_input.Nae = 150;
struct_input.MgATP = 10;
struct_input.MgADP = 0;
struct_input.Pi_total = 0;
struct_input.H_conc = 10^(-4.4);
struct_input.V = -0.08;

Nae_vec = [50; 100; 150];
Vm_vec = struct_Nakao_Gadsby_V_data.V_m;

total_error = 0;
for i_Nae = 1:length(Nae_vec)
    Nae = Nae_vec(i_Nae);
    eval(['cycling_velocity_data_vec = struct_Nakao_Gadsby_V_data.cyc_Nae'...
        num2str(Nae) ';']);
    % Change Nae concentration
    struct_input.Nae = Nae;
    struct_input.V = 0.04;
    normalising_factor = 55/NaK_fitting_vcyc(params_vec,struct_input);
    for i_Vm = 1:length(Vm_vec)
        V_m = Vm_vec(i_Vm);
        % Change membrane voltage
        struct_input.V = V_m;

        weight = 1;
        if Nae >= 100
            weight = weight*4;
        end
        if V_m >= -80/1000 && V_m <= 20/1000
            weight = weight*3;
        end
        
        cycling_velocity_data = cycling_velocity_data_vec(i_Vm);
        cycling_velocity_model = NaK_fitting_vcyc(params_vec,struct_input);
        if Nae ~= 150
            cycling_velocity_model = normalising_factor*cycling_velocity_model;
        end
        
        error = (cycling_velocity_data-cycling_velocity_model)^2;
        total_error = total_error + weight*error;
    end
end

final_error = total_error;