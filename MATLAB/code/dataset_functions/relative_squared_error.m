function answer = relative_squared_error(data_value,model_value)
    relative_error = (model_value-data_value)/data_value;
    answer = relative_error^2;
end