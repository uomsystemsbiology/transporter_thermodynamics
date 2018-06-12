function struct_output = load_Hansen_Nai(file_path)

file_id = fopen(file_path);
data = textscan(file_id,'%s','delimiter','\n');
fclose(file_id);

Nai_vec = [5; 10; 30; 50; 80];
num_Nai_conc = length(Nai_vec);

data = data{1};

struct_output = struct('Nai',Nai_vec,...
    'cyc_rate',zeros(num_Nai_conc,1));

for i_line = 1:num_Nai_conc
    line = data{1+i_line};
    line_split = regexp(line,',','split');
    struct_output.cyc_rate(i_line) = str2double(line_split(2));
end

end