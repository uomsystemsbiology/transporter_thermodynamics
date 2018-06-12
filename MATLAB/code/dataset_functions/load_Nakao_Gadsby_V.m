function struct_output = load_Nakao_Gadsby_V(file_path)

file_id = fopen(file_path);
data = textscan(file_id,'%s','delimiter','\n');
fclose(file_id);

num_Nae_conc = 4;

data = data{1};
num_lines = length(data);
num_values = (num_lines-1)/num_Nae_conc;

struct_output = struct('V_m',transpose((-120:20:60))/1000,...
    'cyc_Nae1p5',zeros(num_values,1),...
    'cyc_Nae50',zeros(num_values,1),...
    'cyc_Nae100',zeros(num_values,1),...
    'cyc_Nae150',zeros(num_values,1));

for i_line = 1:num_values
    line = data{1+num_Nae_conc*i_line};
    line_split = regexp(line,',','split');
    struct_output.cyc_Nae1p5(i_line) = str2double(line_split(2));
    struct_output.cyc_Nae50(i_line) = str2double(line_split(3));
    struct_output.cyc_Nae100(i_line) = str2double(line_split(4));
    struct_output.cyc_Nae150(i_line) = str2double(line_split(5));
end

end