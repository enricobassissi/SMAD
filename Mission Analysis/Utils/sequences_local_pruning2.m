function [asteroid_names, PermutationMatrix_after, HowMany_after] = sequences_local_pruning2(data_elements,p_number)

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
a_asteroids = str2double(data_elements(:,2));
e_asteroids = str2double(data_elements(:,3));
incl_asteroids = str2double(data_elements(:,4));
OM_asteroids = str2double(data_elements(:,5));
om_asteroids = str2double(data_elements(:,6));

HowMany = factorial(length(asteroid_names)) / factorial(length(asteroid_names) - p_number);
[PermutationMatrix_whole, ~] = permnUnique(asteroid_names, p_number);

PermutationMatrix_to_cut = PermutationMatrix_whole;

delta_incl_max = 0.5;
% delta_a_max = 0.5;
for i = 1:HowMany % rows
    idx_first_ast_considered = PermutationMatrix_whole(i,1)==asteroid_names;
    idx_second_ast_considered = PermutationMatrix_whole(i,2)==asteroid_names;
    a_first_ast = a_asteroids(idx_first_ast_considered);
    e_first_ast = e_asteroids(idx_first_ast_considered);
    e_second_ast = e_asteroids(idx_second_ast_considered);
    incl_first_ast = incl_asteroids(idx_first_ast_considered);
    incl_second_ast = incl_asteroids(idx_second_ast_considered);
    if abs(a_first_ast - 1) < 0.2
        if e_first_ast < 0.2
            if incl_first_ast < 2
                sign_of_that_sequence_i = sign(incl_second_ast - incl_first_ast); % if +1 growing i, if -1, decreasing i
                sign_of_that_sequence_e = sign(e_second_ast - e_first_ast); % if +1 growing i, if -1, decreasing i
                for j = 2:p_number % cols
                    idx_ast_considered_a = find(PermutationMatrix_whole(i,j-1)==asteroid_names);
                    idx_ast_considered_b = find(PermutationMatrix_whole(i,j)==asteroid_names);
                    incl_a = incl_asteroids(idx_ast_considered_a);
                    incl_b = incl_asteroids(idx_ast_considered_b);
                    delta_incl = abs(incl_b - incl_a);
                    sign_delta_incl = sign(incl_b - incl_a);
                    if delta_incl > delta_incl_max || sign_delta_incl ~= sign_of_that_sequence_i
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
                    
                    sign_delta_e = sign(ecc_b - ecc_a);
                    % || sign_delta_e ~= sign_of_that_sequence_e
                    if cond_1 || cond_2 
                        PermutationMatrix_to_cut(i,:) = "TO BE CUT w_up";
                    end

            %         a_a = a_asteroids(idx_ast_considered_a);
            %         a_b = a_asteroids(idx_ast_considered_b);
            %         delta_a = abs(a_b - a_a);
            %         if delta_a > delta_a_max
            %             PermutationMatrix_to_cut(i,:) = "TO BE CUT a";
            %         end
                end
            else
                PermutationMatrix_to_cut(i,:) = "TO BE CUT i";
            end
        else
            PermutationMatrix_to_cut(i,:) = "TO BE CUT e";
        end
    else
        PermutationMatrix_to_cut(i,:) = "TO BE CUT a";
    end
end

PermutationMatrix_without_i = PermutationMatrix_to_cut(~all(PermutationMatrix_to_cut == "TO BE CUT i", 2),:); % 2 works by rows, 1 works by columns
PermutationMatrix_without_w_up = PermutationMatrix_without_i(~all(PermutationMatrix_without_i == "TO BE CUT w_up", 2),:); % 2 works by rows, 1 works by columns
PermutationMatrix_without_e = PermutationMatrix_without_w_up(~all(PermutationMatrix_without_w_up == "TO BE CUT e", 2),:); % 2 works by rows, 1 works by columns
PermutationMatrix_after = PermutationMatrix_without_e(~all(PermutationMatrix_without_e == "TO BE CUT a", 2),:); % 2 works by rows, 1 works by columns

HowMany_after = length(PermutationMatrix_after);

end