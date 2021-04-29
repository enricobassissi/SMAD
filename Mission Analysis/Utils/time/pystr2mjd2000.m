function [mjd2000] = pystr2mjd2000(date_with_hyphen)

    spaced_date = replace(date_with_hyphen,"-"," ");
    separated_date = str2double(strsplit(spaced_date," "));

    if numel(separated_date) == 3
        six_elements_vector_date = [separated_date,0,0,0];
    elseif numel(separated_date) == 6
        six_elements_vector_date = separated_date;
    else
        fprintf('insert a date like 2022-01-01 or 2022-01-01-12-35-02');
    end

    mjd2000 = date2mjd2000(six_elements_vector_date);

end

