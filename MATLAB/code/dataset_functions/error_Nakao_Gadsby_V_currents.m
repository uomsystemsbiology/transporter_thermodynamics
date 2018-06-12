function final_error = error_Nakao_Gadsby_V_currents(params_vec,struct_Nakao_Gadsby_I_data)

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
Vm_vec = struct_Nakao_Gadsby_I_data.V_m;

pump_density = params_vec(15);

total_error = 0;
for i_Nae = 1:length(Nae_vec)
    Nae = Nae_vec(i_Nae);
    I_data_vec = struct_Nakao_Gadsby_I_data.I(:,i_Nae);
    % Change Nae concentration
    struct_input.Nae = Nae;
    struct_input.V = 0.04;
    for i_Vm = 1:length(Vm_vec)
        V_m = Vm_vec(i_Vm);
        % Change membrane voltage
        struct_input.V = V_m;

        weight = 1;
        if Nae == 150 && V_m >= 0
            weight = 0;
        end
        
        I_data = I_data_vec(i_Vm);
        I_model = NaK_fitting_vcyc(params_vec,struct_input) * pump_density * 1.6e-5;
        
        error = (I_data-I_model)^2;
        total_error = total_error + weight*error;
    end
end

final_error = total_error;