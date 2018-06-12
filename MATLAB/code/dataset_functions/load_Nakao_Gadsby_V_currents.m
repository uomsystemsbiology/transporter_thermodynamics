function struct_output = load_Nakao_Gadsby_V_currents(file_path)

data = csvread(file_path,1,0);

C_m = 183;

struct_output = struct('V_m',transpose((-140:20:60))/1000,...
    'Nae',[50; 100; 150],...
    'I',data(1:3:end,2:end)/183);

end