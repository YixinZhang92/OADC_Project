function pol_data = convert_data_to_polarity(data)
% This function extracts the polarity of a number, vector or matrix
%
siz_data = size(data);
len_data = siz_data(1)*siz_data(2);

data_flattened = reshape(data,[1 len_data]);

pol_data_tmp = data_flattened./abs(data_flattened);
pol_data = reshape(pol_data_tmp,siz_data);