function an = atan_S_wave_pol(data_x,data_y)
% This function calculates atan of data_x/data_y with angles from 0 to 360
% degrees. The angle is measured clockwisely from the 'north'.
%
% data_x and data_y must be the same size and can be any 2D dimension.
%
siz_data = size(data_x);
len_data = siz_data(1)*siz_data(2);
%
x_flattened = reshape(data_x,[1 len_data]);
y_flattened = reshape(data_y,[1 len_data]);
%
%  zero the an_tmp array
an_tmp=zeros(1,len_data);
%
for i = 1:len_data
    if x_flattened(i) >= 0 && y_flattened(i) >= 0
        x=abs(x_flattened(i)); y= abs(y_flattened(i));
        an_tmp(i) = atand(x/y);
    elseif x_flattened(i) > 0 && y_flattened(i) < 0
        x=abs(x_flattened(i)); y= abs(y_flattened(i));
        an_tmp(i) = 180 - atand(x/y);
    elseif x_flattened(i) <= 0 && y_flattened(i) <= 0
        x=abs(x_flattened(i)); y= abs(y_flattened(i));
        an_tmp(i) = 180 + atand(x/y);
    elseif x_flattened(i) < 0 && y_flattened(i) > 0
        x=abs(x_flattened(i)); y= abs(y_flattened(i));
        an_tmp(i) = 360 - atand(x/y);
    end
end
%
an = reshape(an_tmp,siz_data);