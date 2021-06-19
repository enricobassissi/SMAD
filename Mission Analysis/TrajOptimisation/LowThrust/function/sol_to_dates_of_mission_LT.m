function [sol_dates] = sol_to_dates_of_mission_LT(sol,id_case)

switch id_case
    case 'ga'
        sol_dates.dep = sol.dep_date(1:3);
        date_ga = mjd20002date(sol.MJD0 + sol.TOF0);
        sol_dates.ga = date_ga(1:3);
        date_1ast = mjd20002date(sol.MJD0 + sol.TOF0 + sol.TOF1);
        sol_dates.ast1 = date_1ast(1:3);
        date_2ast = mjd20002date(sol.MJD0 + sol.TOF0 + sol.TOF1 + sol.TOF2);
        sol_dates.ast2 = date_2ast(1:3);
        date_3ast = mjd20002date(sol.MJD0 + sol.TOF0 + sol.TOF1 + sol.TOF2 + sol.TOF3);
        sol_dates.ast3 = date_3ast(1:3);
        date_4ast = mjd20002date(sol.MJD0 + sol.TOF0 + sol.TOF1 + sol.TOF2 + sol.TOF3 + sol.TOF4);
        sol_dates.ast4 = date_4ast(1:3);
        
    case '1RL'
        sol_dates.dep = sol.departure_date(1:3);
        date_1ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOF1);
        sol_dates.ast1_arr = date_1ast_arr(1:3);
        date_1ast_dep = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.CT1);
        sol_dates.ast1_dep = date_1ast_dep(1:3);
        date_2ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.CT1 + sol.TOF2);
        sol_dates.ast2_arr = date_2ast_arr(1:3);
        date_2ast_dep = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.CT1 + ...
            sol.TOF2 + sol.CT2);
        sol_dates.ast2_dep = date_2ast_dep(1:3);
        date_3ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.CT1 + ...
            sol.TOF2 + sol.CT2 + sol.TOF3);
        sol_dates.ast3_arr = date_3ast_arr(1:3);
        date_3ast_dep = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.CT1 + ...
            sol.TOF2 + sol.CT2 + sol.TOF3 + sol.CT3);
        sol_dates.ast3_dep = date_3ast_dep(1:3);
        date_4ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.CT1 + ...
            sol.TOF2 + sol.CT2 + sol.TOF3 + sol.CT3 + sol.TOF4);
        sol_dates.ast4_arr = date_4ast_arr(1:3);
        
    case '1FL'
        sol_dates.dep = sol.departure_date(1:3);
        date_1ast = mjd20002date(sol.MJD0 + sol.TOF1_DIM);
        sol_dates.ast1 = date_1ast(1:3);
        date_2ast = mjd20002date(sol.MJD0 + sol.TOF1_DIM + sol.TOF2_DIM);
        sol_dates.ast2 = date_2ast(1:3);
        date_3ast = mjd20002date(sol.MJD0 + sol.TOF1_DIM + sol.TOF2_DIM + sol.TOF3_DIM);
        sol_dates.ast3 = date_3ast(1:3);
        date_4ast = mjd20002date(sol.MJD0 + sol.TOF1_DIM + sol.TOF2_DIM + sol.TOF3_DIM + sol.TOF4_DIM);
        sol_dates.ast4 = date_4ast(1:3);
        
    case '2RL'
        sol_dates.dep = sol.departure_date(1:3);
        
        date_1ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOF1);
        sol_dates.ast1_arr = date_1ast_arr(1:3);
        date_1ast_dep = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.CT1);
        sol_dates.ast1_dep = date_1ast_dep(1:3);
        date_2ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.CT1 + sol.TOF2);
        sol_dates.ast2_arr = date_2ast_arr(1:3);
        
        date_a_ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOFa);
        sol_dates.asta_arr = date_a_ast_arr(1:3);
        date_a_ast_dep = mjd20002date(sol.departure_mjd2000 + sol.TOFa + sol.CTa);
        sol_dates.asta_dep = date_a_ast_dep(1:3);
        date_b_ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOFa + sol.CTa + ...
            sol.TOFb);
        sol_dates.astb_arr = date_b_ast_arr(1:3);
	case '2RL_GA_in_between_asteroids'
        sol_dates.dep = sol.departure_date(1:3);
        
        date_1ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOF1);
        sol_dates.ast1_arr = date_1ast_arr(1:3);
        date_1ast_dep = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.CT1);
        sol_dates.ast1_dep = date_1ast_dep(1:3);
        date_GA1 = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.CT1 + sol.TOFGA1);
        sol_dates.GA1 = date_GA1(1:3);
        date_2ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.CT1 + sol.TOFGA1 + sol.TOF2);
        sol_dates.ast2_arr = date_2ast_arr(1:3);
        
        date_a_ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOFa);
        sol_dates.asta_arr = date_a_ast_arr(1:3);
        date_a_ast_dep = mjd20002date(sol.departure_mjd2000 + sol.TOFa + sol.CTa);
        sol_dates.asta_dep = date_a_ast_dep(1:3);
        date_GAa = mjd20002date(sol.departure_mjd2000 + sol.TOFa + sol.CTa + sol.TOFGAa);
        sol_dates.GAa = date_GAa(1:3);
        date_b_ast_arr = mjd20002date(sol.departure_mjd2000 + sol.TOFa + sol.CTa + ...
            sol.TOFGAa + sol.TOFb);
        sol_dates.astb_arr = date_b_ast_arr(1:3);
        
    case '2FL'
        sol_dates.dep = sol.departure_date(1:3);
        
        date_1ast = mjd20002date(sol.departure_mjd2000 + sol.TOF1);
        sol_dates.ast1 = date_1ast(1:3);
        date_2ast = mjd20002date(sol.departure_mjd2000 + sol.TOF1 + sol.TOF2);
        sol_dates.ast2 = date_2ast(1:3);
        
        date_aast = mjd20002date(sol.departure_mjd2000 + sol.TOFa);
        sol_dates.asta = date_aast(1:3);
        date_bast = mjd20002date(sol.departure_mjd2000 + sol.TOFa + sol.TOFb);
        sol_dates.astb = date_bast(1:3);
        
end

end