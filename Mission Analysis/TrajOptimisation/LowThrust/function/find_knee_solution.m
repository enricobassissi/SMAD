function [knee_idx, d] = find_knee_solution(Fval)

    d = zeros(size(Fval,1),1);
    max_obj1 = max(Fval(:,1));
    max_obj2 = max(Fval(:,2));
    for i = 1:size(Fval,1)
        d(i,1) = sqrt((Fval(i,1)/max_obj1)^2+(Fval(i,2)/max_obj2)^2);
    end

    knee_idx = find(min(d)==d);

end