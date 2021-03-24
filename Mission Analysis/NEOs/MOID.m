function [min_moid, date_min_moid] = MOID(r, r_EA, date_vect)

%{
OUTPUT:
    min_moid: AU
    date_min_moid: date_format 6x1
%}
AU = astroConstants(2);

moid = zeros(length(r),1);
for i = 1:length(r)
    moid(i,1) = abs(norm(r(i,:)) - norm(r_EA(i,:)));
end

min_moid = min(moid)/AU;

[idx_min_moid] = find(min(moid)==moid);
date_min_moid = date_vect(idx_min_moid,:);

end