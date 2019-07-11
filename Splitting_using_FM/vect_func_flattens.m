function m = vect_func_flattens(m)

siz_data = size(m);

if length(siz_data) == 3;
    len_data = siz_data(1)*siz_data(2)*siz_data(3);
elseif length(siz_data) == 2;
    len_data = siz_data(1)*siz_data(2);
end

m = reshape(m,[1 len_data]);