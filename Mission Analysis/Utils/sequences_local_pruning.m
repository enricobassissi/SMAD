function [asteroid_names, PermutationMatrix_after, HowMany_after] = sequences_local_pruning(data_elements,p_number)

%{
 p_number = number of selected element in the dataset
 w_up = Longitude of pericenter, OM+om
 transfer trajectory between two orbits (subscripts 1 and 2) with a large 
 delta_w_up = mod((w_up_1 - w_up_2), pi) is, in fact, as more difficult to 
 achieve as the eccentricities of the two orbits increase. 
 For this reason, a threshold on the maximum variation of the longitude of 
 pericenter has been considered for each object, taking into account the 
 value of the eccentricity as follows: delta_w_up_max := pi(1-e)^2
%}

% asteroids elements extraction
asteroid_names = data_elements(:,1);
e_asteroids = str2double(data_elements(:,3));
incl_asteroids = str2double(data_elements(:,4));
OM_asteroids = str2double(data_elements(:,5));
om_asteroids = str2double(data_elements(:,6));

HowMany = factorial(length(asteroid_names)) / factorial(length(asteroid_names) - p_number);
[PermutationMatrix_whole, ~] = permnUnique(asteroid_names, p_number);

PermutationMatrix_to_cut = PermutationMatrix_whole;

delta_incl_max = 5;
for i = 1:HowMany % rows
    for j = 2:p_number % cols
        idx_ast_considered_a = find(PermutationMatrix_whole(i,j-1)==asteroid_names);
        idx_ast_considered_b = find(PermutationMatrix_whole(i,j)==asteroid_names);
        incl_a = incl_asteroids(idx_ast_considered_a);
        incl_b = incl_asteroids(idx_ast_considered_b);
        delta_incl = abs(incl_b - incl_a);
        if delta_incl > delta_incl_max
            PermutationMatrix_to_cut(i,:) = "TO BE CUT i";
        end
        
        w_up_a = deg2rad(OM_asteroids(idx_ast_considered_a) + om_asteroids(idx_ast_considered_a));
        w_up_b = deg2rad(OM_asteroids(idx_ast_considered_b) + om_asteroids(idx_ast_considered_b));
%         delta_w_up = w_up_b - w_up_a;
        ecc_a = e_asteroids(idx_ast_considered_a);
        ecc_b = e_asteroids(idx_ast_considered_b);
        delta_w_up_max_a = pi*(1-ecc_a)^2; %pi*(1-ecc_a)^2
        delta_w_up_max_b = pi*(1-ecc_b)^2;
        % By using the threshold, the arrival object is removed from the
        % locally pruned database if at least one of the following conditions 
        % is satisfied:
        cond_1 = mod(w_up_a+delta_w_up_max_a+pi,2*pi) < mod(w_up_b-delta_w_up_max_b+pi,2*pi);
        cond_2 = mod(w_up_b+delta_w_up_max_b+pi,2*pi) < mod(w_up_a-delta_w_up_max_a+pi,2*pi);
        if cond_1 || cond_2
            PermutationMatrix_to_cut(i,:) = "TO BE CUT w_up";
        end
%         if (incl_asteroids(idx_ast_considered_b) - incl_asteroids(idx_ast_considered_a)) > 5 && (e_asteroids(idx_ast_considered_b) - e_asteroids(idx_ast_considered_a)) > 0.2
%             PermutationMatrix(i,:) = "TO BE CUT i and e";
%         end
    end
end

PermutationMatrix_without_i = PermutationMatrix_to_cut(~all(PermutationMatrix_to_cut == "TO BE CUT i", 2),:); % 2 works by rows, 1 works by columns
PermutationMatrix_after = PermutationMatrix_without_i(~all(PermutationMatrix_without_i == "TO BE CUT w_up", 2),:); % 2 works by rows, 1 works by columns

HowMany_after = length(PermutationMatrix_after);

end