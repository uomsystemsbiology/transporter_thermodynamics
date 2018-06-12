function struct_output = load_Nakao_Gadsby_Ke(file_path)

file_id = fopen(file_path);
data = textscan(file_id,'%s','delimiter','\n');
fclose(file_id);

Ke_vec = [0.1; 0.3; 1.0; 2.7; 5.4; 10];
num_Ke_conc = length(Ke_vec);

data = data{1};

struct_output = struct('Ke',Ke_vec,...
    'cyc_rate',zeros(num_Ke_conc,1));

for i_line = 1:num_Ke_conc
    line = data{1+i_line};
    line_split = regexp(line,',','split');
    struct_output.cyc_rate(i_line) = str2double(line_split(2));
end

% Set normalised cycling rate to 1 for [Ke] = 5.4mM
struct_output.cyc_rate(5) = 1;

end