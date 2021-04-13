function [ref_sel_orb_el, units] = get_orbital_elements_and_sigma(asteroid_list,py_refined_selected_dict)

%dict_selected_asteroids = cell(length(asteroid_list),1);
ref_sel_orb_el = cell(length(asteroid_list),1);

for name = 1:length(asteroid_list)
    % py dict extraction
    py_dict_selected_asteroids = py_refined_selected_dict{asteroid_list(name)};
    %dict_selected_asteroids{name} = py_dict_selected_asteroids; 
    %struct_sel_ast(name) = struct(struct(dict_selected_asteroids{name})); %= cellfun(@string,cell(py_single_ast_dict));
    
    % Extraction of orbital param of the orbits with uncertainty
    a = py_dict_selected_asteroids{'orbital_elements'}{'a'};
    a_sig = py_dict_selected_asteroids{'orbital_elements'}{'a_sig'};
    
    e = py_dict_selected_asteroids{'orbital_elements'}{'e'};
    e_sig = py_dict_selected_asteroids{'orbital_elements'}{'e_sig'};
    
    i = py_dict_selected_asteroids{'orbital_elements'}{'i'};
    i_sig = py_dict_selected_asteroids{'orbital_elements'}{'i_sig'};
    
    OM = py_dict_selected_asteroids{'orbital_elements'}{'om'};
    OM_sig = py_dict_selected_asteroids{'orbital_elements'}{'om_sig'};
    
    om = py_dict_selected_asteroids{'orbital_elements'}{'w'};
    om_sig = py_dict_selected_asteroids{'orbital_elements'}{'w_sig'};
    
    M = py_dict_selected_asteroids{'orbital_elements'}{'ma'};
    M_sig = py_dict_selected_asteroids{'orbital_elements'}{'ma_sig'};
    
    ref_sel_orb_el{name} = [a.scale,a_sig.scale; 
                            str2double(string(e)), str2double(string(e_sig)); 
                            i.scale,i_sig.scale; 
                            OM.scale,OM_sig.scale; 
                            om.scale,om_sig.scale; 
                            M.scale,M_sig.scale];
end

units = {'AU'; '-'; 'deg'; 'deg'; 'deg'; 'deg'};

end