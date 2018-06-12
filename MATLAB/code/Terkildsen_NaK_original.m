function v_cyc_nak = Terkildsen_NaK_original(struct_input)
% Calculate the cycling velocity of the Na/K ATPase using parameters and
% expressions from the original model.

F = 96485; % unit C/mol
R = 8.314; % unit J/mol/K
T = struct_input.T; % unit K

Vm = struct_input.V; % Unit V

% Concentrations in mM
MgATP = struct_input.MgATP;
Pi_total = struct_input.Pi_total;
Nai = struct_input.Nai;
Nae = struct_input.Nae;
Ke = struct_input.Ke;
Ki = struct_input.Ki;
H_conc = struct_input.H_conc;
MgADP = struct_input.MgADP;

K0d_nak_nai= 0.1209868E-07; 		%K0d_nak_nai (mM)
K0d_nak_nae= 141.271240234375000; 	%K0d_nak_nae (mM)
Kd_nak_nai2= 154.389923095703125; 	%Kd_nak_nai2 (mM)
Kd_nak_nae2= 356.009094238281250 ;	%Kd_nak_nae2 (mM)
Kd_nak_ki= 255.130677771665574 ;	%Kd_nak_ki (mM)
Kd_nak_ke= 6.80234527587890625 ;	%Kd_nak_ke (mM)
Kd_nak_MgATP= 2.92630076408386230 ;	%K_nak_mgatp (mM)
factor_nak_del= -0.179529249668121338; 	%factor_nak_del
rate_f_nak_1= 1664.12402343750000; 	%rate_f_1 (s-1)
rate_f_nak_2= 110.425056457519531 ;	%rate_f_2 (s-1)
rate_f_nak_3= 2314.25756835937500 ;	%rate_f_3 (s-1)
rate_f_nak_4= 462.377563476562500 ;	%rate_f_4 (s-1)
rate_r_nak_1= 264.146662389765481 ;	%rate_r_1 (s-1 mM-1)
rate_r_nak_2= 14.0371952056884766 ;	%rate_r_2  (s-1 mM-1)
rate_r_nak_3= 7599999.50000000000 ;	%rate_r_3 (s-1)
rate_r_nak_4= 1702.52624511718750 ;	%rate_r_4  (s-1)

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
 		(Nai_bound2^2.0d0)/(((1.0d0+Nai_bound1)*...
 		(1.0d0+Nai_bound2)^2.0d0)+...
            ((1.0d0+Ki_bound)^2.0d0)-1.0d0);

alpha_f_nak_2 = rate_f_nak_2;

alpha_f_nak_3 = rate_f_nak_3*(Ke_bound^2.0d0)/...
 		(((1.0d0+Nae_bound1)*(1.0d0+Nae_bound2)^2.0d0)+...
 		((1.0d0+Ke_bound)^2.0d0)-1.0d0);

alpha_f_nak_4 = rate_f_nak_4*MgATP_bound/...
 		(1.0d0+MgATP_bound);

%Reverse rate constants
alpha_r_nak_1 = rate_r_nak_1*MgADP;

alpha_r_nak_2 = rate_r_nak_2*Nae_bound1*...
		((Nae_bound2)^2.0d0)/...
 		(((1.0d0+Nae_bound1)*(1.0d0+Nae_bound2)^2.0d0)+...
		((1.0d0+Ke_bound)^2.0d0)-1.0d0);

alpha_r_nak_3 = rate_r_nak_3*Pi*H_conc/...
 		(1.0d0+MgATP_bound);

alpha_r_nak_4 = rate_r_nak_4*(Ki_bound^2.0d0)/...
 		(((1.0d0+Nai_bound1)*(1.0d0+Nai_bound2)^2.0d0)+...
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