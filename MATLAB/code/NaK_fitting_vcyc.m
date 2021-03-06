function v_cyc_nak = NaK_fitting_vcyc(params,struct_input)
% Calculate cycling velocity of the updated model from parameters

F = 96485; % unit C/mol
R = 8.314; % unit J/mol/K
T = struct_input.T; % unit K

G_MgATP_0 = 11900;
T_body = 310;
K_MgATP_hyd = exp(-G_MgATP_0/R/T_body)*10^6;

Vm = struct_input.V; % Unit volts

% Concentrations in mM
MgATP = struct_input.MgATP;
Pi_total = struct_input.Pi_total;
Nai = struct_input.Nai;
Nae = struct_input.Nae;
Ke = struct_input.Ke;
Ki = struct_input.Ki;
H_conc = struct_input.H_conc;
MgADP = struct_input.MgADP;

K0d_nak_nai= params(1); 		%K0d_nak_nai (mM)
K0d_nak_nae= params(2); 	%K0d_nak_nae (mM)
Kd_nak_nai2= params(3); 	%Kd_nak_nai2 (mM)
Kd_nak_nae2= params(4);	%Kd_nak_nae2 (mM)
Kd_nak_ke= params(5);	%Kd_nak_ke (mM)
Kd_nak_MgATP= params(6);	%K_nak_mgatp (mM)
rate_f_nak_1= params(7); 	%rate_f_1 (s-1)
rate_f_nak_2= params(8);	%rate_f_2 (s-1)
rate_f_nak_3= params(9);	%rate_f_3 (s-1)
rate_f_nak_4= params(10);	%rate_f_4 (s-1)
rate_r_nak_1= rate_f_nak_1/6.3;	%rate_r_1 (s-1 mM-1)
rate_r_nak_2= params(11);	%rate_r_2  (s-1 mM-1)
rate_r_nak_3= params(12);	%rate_r_3 (s-1)
rate_r_nak_4= params(13);	%rate_r_4  (s-1)
Kd_nak_ki= sqrt(rate_r_nak_1*rate_r_nak_2*rate_r_nak_3*rate_r_nak_4/...
    rate_f_nak_1/rate_f_nak_2/rate_f_nak_3/rate_f_nak_4*...
    K0d_nak_nai*Kd_nak_nai2^2*Kd_nak_ke^2*Kd_nak_MgATP*K_MgATP_hyd/...
    K0d_nak_nae/Kd_nak_nae2^2);	%Kd_nak_ki (mM)
factor_nak_del= params(14); 	%factor_nak_del

pKd_nak_hpi = 6.77d0;	%pKd_nak_hpi
Kd_nak_kpi = 292.0d0;	%Kd_nak_kpi (mM)
Kd_nak_napi = 224.0d0;	%Kd_nak_napi (mM)

Kd_nak_nae1 = K0d_nak_nae*exp((1.0d0 ...
    +factor_nak_del)*F*Vm/(R*T));

Kd_nak_nai1 = K0d_nak_nai*exp(factor_nak_del*...
    F*Vm/(R*T));

Nae_bound1 = Nae/Kd_nak_nae1;

Nai_bound1 = Nai/Kd_nak_nai1;

Nae_bound2 = Nae/Kd_nak_nae2;

Nai_bound2 = Nai/Kd_nak_nai2;

Ke_bound = Ke/Kd_nak_ke;

Ki_bound = Ki/Kd_nak_ki;

MgATP_bound = MgATP/Kd_nak_MgATP;

Pi = Pi_total/(1.0d0+(Ki/Kd_nak_kpi)+ ...
 	(H_conc/(10^(3.0-pKd_nak_hpi)))+(Nai/...
 	Kd_nak_napi));


%Forward rate constants
alpha_f_nak_1 = rate_f_nak_1*Nai_bound1*...
 		(Nai_bound2^2.0d0)/(Nai_bound1*Nai_bound2^2 + ...
        ((1.0d0+Nai_bound2)^2.0d0)+...
            ((1.0d0+Ki_bound)^2.0d0)-1.0d0);

alpha_f_nak_2 = rate_f_nak_2;

alpha_f_nak_3 = rate_f_nak_3*(Ke_bound^2.0d0)/...
 		(Nae_bound1*Nae_bound2^2+...
        ((1.0d0+Nae_bound2)^2.0d0)+...
 		((1.0d0+Ke_bound)^2.0d0)-1.0d0);

alpha_f_nak_4 = rate_f_nak_4*MgATP_bound/...
 		(1.0d0+MgATP_bound);

%Reverse rate constants
alpha_r_nak_1 = rate_r_nak_1*MgADP;

alpha_r_nak_2 = rate_r_nak_2*Nae_bound1*...
		((Nae_bound2)^2.0d0)/...
 		(Nae_bound1*Nae_bound2^2 +...
        (1.0d0+Nae_bound2)^2.0d0+...
		((1.0d0+Ke_bound)^2.0d0)-1.0d0);

alpha_r_nak_3 = rate_r_nak_3*Pi*H_conc/...
 		(1.0d0+MgATP_bound);

alpha_r_nak_4 = rate_r_nak_4*(Ki_bound^2.0d0)/...
 		(Nai_bound1*Nai_bound2^2+...
        ((1.0d0+Nai_bound2)^2.0d0)+...
		((1.0d0+Ki_bound)^2.0d0)-1.0d0);


sig_nak = alpha_r_nak_1*alpha_r_nak_2*alpha_r_nak_3 + ...
 		alpha_r_nak_1*alpha_r_nak_2*alpha_f_nak_4 + ...
 		alpha_r_nak_1*alpha_f_nak_3*alpha_f_nak_4 + ...
 		alpha_f_nak_2*alpha_f_nak_3*alpha_f_nak_4 + ...
 		alpha_r_nak_2*alpha_r_nak_3*alpha_r_nak_4 + ...
 		alpha_f_nak_1*alpha_r_nak_2*alpha_r_nak_3 +  ...
 		alpha_f_nak_1*alpha_r_nak_2*alpha_f_nak_4 +  ...
 		alpha_f_nak_1*alpha_f_nak_3*alpha_f_nak_4 + ...
 		alpha_r_nak_1*alpha_r_nak_3*alpha_r_nak_4 + ...
 		alpha_f_nak_2*alpha_r_nak_3*alpha_r_nak_4 + ...
 		alpha_f_nak_1*alpha_f_nak_2*alpha_r_nak_3 +  ...
 		alpha_f_nak_1*alpha_f_nak_2*alpha_f_nak_4 +  ...
 		alpha_r_nak_1*alpha_r_nak_2*alpha_r_nak_4 +  ...
 		alpha_r_nak_1*alpha_f_nak_3*alpha_r_nak_4 + ...
 		alpha_f_nak_2*alpha_f_nak_3*alpha_r_nak_4 + ...
 		alpha_f_nak_1*alpha_f_nak_2*alpha_f_nak_3;

v_cyc_nak = (alpha_f_nak_1*alpha_f_nak_2*alpha_f_nak_3*...
 		alpha_f_nak_4-alpha_r_nak_1*alpha_r_nak_2*...
 		alpha_r_nak_3*alpha_r_nak_4)/sig_nak;

end