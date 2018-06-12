function total_error = NaK_total_error(params_vec,struct_NaK_data)

struct_Nakao_Gadsby_V_data = struct_NaK_data.Nakao_Gadsby_V_data;
struct_Nakao_Gadsby_I_data = struct_NaK_data.Nakao_Gadsby_I_data;
struct_Nakao_Gadsby_Ke_data = struct_NaK_data.Nakao_Gadsby_Ke_data;
struct_Hansen_Nai_data = struct_NaK_data.Hansen_Nai_data;
struct_Friedrich_MgATP_data = struct_NaK_data.Friedrich_MgATP_data;

weight_Nakao_Gadsby_V = 1;
weight_Nakao_Gadsby_V_currents = 1000;
weight_Nakao_Gadsby_Ke = 30000;
weight_Hansen_Nai = 100000;
weight_Friedrich_MgATP = 100000;

total_error = weight_Nakao_Gadsby_V*error_Nakao_Gadsby_V_unc(params_vec,struct_Nakao_Gadsby_V_data) ...
    + weight_Nakao_Gadsby_V_currents*error_Nakao_Gadsby_V_currents(params_vec,struct_Nakao_Gadsby_I_data) ...
    + weight_Nakao_Gadsby_Ke*error_Nakao_Gadsby_Ke(params_vec,struct_Nakao_Gadsby_Ke_data) ...
    + weight_Hansen_Nai*error_Hansen_Nai(params_vec,struct_Hansen_Nai_data) ...
    + weight_Friedrich_MgATP*error_Friedrich_ATP(params_vec,struct_Friedrich_MgATP_data);

end