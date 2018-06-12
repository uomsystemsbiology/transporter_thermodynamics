function struct_output = load_Friedrich_ATP(file_path)
% Loads ATP data from Friedrich et al.

file_id = fopen(file_path);
data = textscan(file_id,'%s','delimiter','\n');
fclose(file_id);

data = data{1};
num_data_points = length(data)-1;

struct_output = struct('MgATP',zeros(num_data_points,1),...
    'cyc_rate',zeros(num_data_points,1));

for i_line = 1:num_data_points
    line = data{1+i_line};
    line_split = regexp(line,',','split');
    struct_output.MgATP(i_line) = str2double(line_split(1))*1e-3;
    struct_output.cyc_rate(i_line) = str2double(line_split(2));
end

end