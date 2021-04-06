from pprint import pprint
import neo_api_function as neo
def main_fun():
    esa_risk_names=neo.extract_esa_name_from_file("esa_risk_list.txt")
    sentry_risk_names=neo.get_sentry_risk_list()
    risk_list=neo.merge_risk_lists(esa_risk_names, sentry_risk_names)
    dict_risk_list=neo.get_dict(risk_list)

    # MOID<=0.05au, H<=26 (if H is not available diameter>=200m)
    MOID_H_selected=[]
    for key in dict_risk_list.keys():
        if float(dict_risk_list[key]['moid'].scale)<=0.05: #MOID<=0.05 AU
            if (dict_risk_list[key]["magn_radius_flag"]=='H' and float(dict_risk_list[key]["H"])<=26) or (dict_risk_list[key]["magn_radius_flag"]=='D' and float(dict_risk_list[key]["D"])>=200):
                MOID_H_selected.append(key)

    # At least one impact 2026<year<2048 with a Palermo Scale>=-7
    date_selected=[]
    PS_date_selected=[]
    for key in dict_risk_list.keys():
        if '0' in dict_risk_list[key]['impacts'].keys():
            max_P=-100;
            date_flag=0;
            for imp_id in dict_risk_list[key]['impacts'].keys():
                word='';
                for c in dict_risk_list[key]['impacts'][imp_id]['date']:
                    #pprint(c)
                    if c=='-':
                        break;
                    else:
                        word=word+c;
                if int(word)<2048 and int(word)>2026:
                    date_flag=1;
                    if float(dict_risk_list[key]['impacts'][imp_id]['ps'])>max_P:
                        max_P=float(dict_risk_list[key]['impacts'][imp_id]['ps']);
            if date_flag==1:
                date_selected.append(key)
                if max_P>=-7:
                    PS_date_selected.append(key)
            dict_risk_list[key]["PS"]=max_P;        

    # Orbit Uncertantains filter (number of observation>=40)
    OU_selected=[]
    for key in dict_risk_list:
            if int(dict_risk_list[key]['N_obs'])>=40:
                OU_selected.append(key)
    #Intersection of filtered lists
    refined_selected=list(set(list(set(PS_date_selected) & set(MOID_H_selected))) & set(OU_selected))
    i=0
    refined_dict={}
    PS_list=[]
    for selected in refined_selected:            
                PS_list.append(dict_risk_list[selected]['PS'])

    index_list=[i[0] for i in sorted(enumerate(PS_list), key=lambda x:x[1])]
    index_list.reverse()
    for ind in index_list:
        pprint(refined_selected[ind])
        pprint('Max Palermo Scale:' + str(PS_list[ind]))
        pprint('OCC:' + str(dict_risk_list[refined_selected[ind]]['condition_code']))
    return(index_list)