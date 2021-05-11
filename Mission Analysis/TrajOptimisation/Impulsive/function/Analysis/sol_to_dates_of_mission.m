function [sol_dates] = sol_to_dates_of_mission(sol,id_case)

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
    case '1SC_RV'
        sol_dates.dep = sol.dep_date(1:3);
        date_1ast_arr = mjd20002date(sol.MJD0 + sol.TOF1);
        sol_dates.ast1_arr = date_1ast_arr(1:3);
        date_1ast_dep = mjd20002date(sol.MJD0 + sol.TOF1 + sol.buffer_time1);
        sol_dates.ast1_dep = date_1ast_dep(1:3);
        date_2ast_arr = mjd20002date(sol.MJD0 + sol.TOF1 + sol.buffer_time1 + sol.TOF2);
        sol_dates.ast2_arr = date_2ast_arr(1:3);
        date_2ast_dep = mjd20002date(sol.MJD0 + sol.TOF1 + sol.buffer_time1 + ...
            sol.TOF2 + sol.buffer_time2);
        sol_dates.ast2_dep = date_2ast_dep(1:3);
        date_3ast_arr = mjd20002date(sol.MJD0 + sol.TOF1 + sol.buffer_time1 + ...
            sol.TOF2 + sol.buffer_time2 + sol.TOF3);
        sol_dates.ast3_arr = date_3ast_arr(1:3);
        date_3ast_dep = mjd20002date(sol.MJD0 + sol.TOF1 + sol.buffer_time1 + ...
            sol.TOF2 + sol.buffer_time2 + sol.TOF3 + sol.buffer_time3);
        sol_dates.ast3_dep = date_3ast_dep(1:3);
        date_4ast_arr = mjd20002date(sol.MJD0 + sol.TOF1 + sol.buffer_time1 + ...
            sol.TOF2 + sol.buffer_time2 + sol.TOF3 + sol.buffer_time3 + sol.TOF4);
        sol_dates.ast4_arr = date_4ast_arr(1:3);
    case '1SC_FB'
        sol_dates.dep = sol.dep_date(1:3);
        date_1ast = mjd20002date(sol.MJD0 + sol.TOF1);
        sol_dates.ast1 = date_1ast(1:3);
        date_2ast = mjd20002date(sol.MJD0 + sol.TOF1 + sol.TOF2);
        sol_dates.ast2 = date_2ast(1:3);
        date_3ast = mjd20002date(sol.MJD0 + sol.TOF1 + sol.TOF2 + sol.TOF3);
        sol_dates.ast3 = date_3ast(1:3);
        date_4ast = mjd20002date(sol.MJD0 + sol.TOF1 + sol.TOF2 + sol.TOF3 + sol.TOF4);
        sol_dates.ast4 = date_4ast(1:3);
    case '2SC_RV'
        sol_dates.dep = sol.dep_date(1:3);
        
        date_1ast_arr = mjd20002date(sol.MJD0 + sol.TOF1);
        sol_dates.ast1_arr = date_1ast_arr(1:3);
        date_1ast_dep = mjd20002date(sol.MJD0 + sol.TOF1 + sol.buffer_time1);
        sol_dates.ast1_dep = date_1ast_dep(1:3);
        date_2ast_arr = mjd20002date(sol.MJD0 + sol.TOF1 + sol.buffer_time1 + sol.TOF2);
        sol_dates.ast2_arr = date_2ast_arr(1:3);
        
        date_a_ast_arr = mjd20002date(sol.MJD0 + sol.TOFa);
        sol_dates.asta_arr = date_a_ast_arr(1:3);
        date_a_ast_dep = mjd20002date(sol.MJD0 + sol.TOFa + sol.buffer_time2);
        sol_dates.asta_dep = date_a_ast_dep(1:3);
        date_b_ast_arr = mjd20002date(sol.MJD0 + sol.TOFa + sol.buffer_time2 + ...
            sol.TOFb);
        sol_dates.astb_arr = date_b_ast_arr(1:3);
    case '2SC_FB'
        sol_dates.dep = sol.dep_date(1:3);
        
        date_1ast = mjd20002date(sol.MJD0 + sol.TOF1);
        sol_dates.ast1 = date_1ast(1:3);
        date_2ast = mjd20002date(sol.MJD0 + sol.TOF1 + sol.TOF2);
        sol_dates.ast2 = date_2ast(1:3);
        
        date_aast = mjd20002date(sol.MJD0 + sol.TOFa);
        sol_dates.asta = date_aast(1:3);
        date_bast = mjd20002date(sol.MJD0 + sol.TOFa + sol.TOFb);
        sol_dates.astb = date_bast(1:3);
        
        
end

end