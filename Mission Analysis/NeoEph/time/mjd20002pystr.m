function [pystr_date] = mjd20002pystr(mjd)

    date_long_six_element_vector = mjd20002date(mjd);

    YMD = string(date_long_six_element_vector(1:3));
    joined_YMD = join(YMD);
    pystr_date = replace(joined_YMD," ","-");

end

