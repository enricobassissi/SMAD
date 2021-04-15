import requests
import math
import numpy as np
from astroquery.jplhorizons import Horizons
from astroquery.jplsbdb import SBDB
from astropy.time import Time
from astropy.table import QTable
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from datetime import datetime
from array import *
import re
import pandas as pd

"""
# This API provides access to the JPL/SSD small-body mission design suite. The following operation modes are available:
# Mode A (accessible) - retrieves the list of accessible small bodies based on user-defined constraint.
# Mode Q (query) - retrieves pre-computed mission options to a specific object. Both impulsive and low-thrust gravity-assist mission options are available.
# Mode M (map) - an extension of mode Q for the impulsive case, returns the data required to render a porkchop plot with multiple parameters.
# Mode T (mission extension) - retrieves the list of small bodies that come closest (or within a prescribed distance) to a user-defined orbit during a certain period of time. This is a crude filter for finding potential candidates for mission extensions.
# This script emphazise the development of the following mode
# Mode M - In addition to querying the database like in mode Q (ballistic), compute all ballistic mission options to the specified target within certain ranges of launch dates and times of flight.
# In addition, the values of the x-y axes of the maps are also provided:
# dep_date - departure dates from Earth (Modified Julian Date), corresponding to the x-axis.
# tof - times of flight to the target (days), corresponding to the y-axis.
# If dep_date has m elements and tof has n, then the 2D arrays are of dimension n x m.
# vinf_dep
# vinf_arr
# Reference of the API construction, data input and output: https://ssd-api.jpl.nasa.gov/doc/mdesign.html
"""

def get_mission_profiles(asteroid_name,mjd0,span,tof_min,tof_max,step):
    """
    # asteroid_name:    designation (provisional or IAU-number) of the desired object (e.g., 2015 AB or 141P or 433).
    #                   NOTE: when submitting a des containing a space in your query string, you must replace the space with %20, for example 2015%20AB
    # mjd0:             first launch date to be explored (Modified Julian Date)
    # span:             duration of the launch-date period to be explored (days)
    # tof-min:          minimum time of flight to be considered (days)
    # tof-max:          maximum time of flight to be considered (days)
    # step:             time step used to advance both the launch date and the time of flight (days). 
    """

    # The size of the transfer map is limited to 1,500,000 points
    sim_lim_points = 1500000 #1.5 millions
    if int(span)/int(step) > sim_lim_points:
        print('outside of tool limits') # TODO return error
    
    # Construction of the HTTP request
    url_base = 'https://ssd-api.jpl.nasa.gov/mdesign.api'
    url = f'{url_base}?sstr={str(asteroid_name)}&mjd0={str(mjd0)}&span={str(span)}&tof-min={str(tof_min)}&tof-max={str(tof_max)}&step={str(step)}'
    r = requests.get(url)
    data = r.json()

    # Elaboration of data
    available_missions = len(data['selectedMissions'])
    mission_profiles={}
    mjd01Jan2021 = 59215
    mjd31Dec2028 = 62136
    mjd01Jan2048 = 69076
    j = 0
    for mission_id in range(available_missions):
        if (data["selectedMissions"][mission_id][0] > mjd01Jan2021 and data["selectedMissions"][mission_id][0] < mjd31Dec2028 and data["selectedMissions"][mission_id][1] < mjd01Jan2048):
            sel_profile={"fullname": data["object"]["fullname"],
                      "mjd0": data["selectedMissions"][mission_id][0],
                      "mjdf": data["selectedMissions"][mission_id][1],
                      "tof": data["selectedMissions"][mission_id][9],
                      "vinf_dep": data["selectedMissions"][mission_id][2],
                      "vinf_arr": data["selectedMissions"][mission_id][3],
                      "earth_dist": data["selectedMissions"][mission_id][5],
                      "phase_ang": data["selectedMissions"][mission_id][4], 
                      "elong_arr": data["selectedMissions"][mission_id][6], 
                      "decl_dep": data["selectedMissions"][mission_id][7],
                      "approach_ang": data["selectedMissions"][mission_id][8],      
                     }
            mission_profiles[str(j)]=sel_profile
            j = j+1
    
    # Find min dv mission profile
    mission_profile_min_dv  = get_min_dv_mission_profile(mission_profiles) #mp_dv_plot removed
    
    # Porkchop data
    porkchop_dv, dep_date, tof = get_mission_porkchop(data) #pc_plot removed
    #return mission_profiles    
    return mission_profiles, porkchop_dv, dep_date, tof,  mission_profile_min_dv #pc_plot, , mp_dv_plot removed

def get_min_dv_mission_profile(mission_profiles):

     # Find the best Mission Profile
    dv = np.zeros(len(mission_profiles.keys()))
    for profile in range(len(mission_profiles)):
        dv[profile] = mission_profiles[str(profile)]['vinf_dep'] + mission_profiles[str(profile)]['vinf_arr']

    index = np.linspace(0,len(dv)-1,len(dv))
    mask = [dv == np.min(dv)]
    idx_min = index[tuple(mask)]
    if type(idx_min) == float:
        mission_profile_min_dv = mission_profiles[str(int(idx_min))]
    else:
        #val, idx_min = min((val, idx) for (idx, val) in enumerate(dv))
        mission_profile_min_dv = mission_profiles[str(int(idx_min[0]))] # NOTE: there could be more than one best solution, here i took arbitrarly the first
    
    # Plot of the mission profiles and highlight the best one
    # fig = plt.figure()
    # plt.plot(dv, "*");
    # fig.suptitle('Mission Profile dv Distribution for '+ mission_profile_min_dv["fullname"])
    # plt.xlabel('$idx$ (-)')
    # plt.ylabel('$dv$ (km/s)')
    # plt.plot(int(idx_min[0]), dv[int(idx_min[0])], "r+");
    
    return mission_profile_min_dv#, fig

def get_mission_porkchop(data):

    # to pass mjd2000 to date format
    # dep_date=Time(data["dep_date"], format='mjd').to_value('iso', 'date');
    dep_date=data["dep_date"]
    tof=data["tof"]
    
    # Elaboration of data
    m = len(dep_date)
    n = len(tof)
    porkchop_map = np.zeros([n,m]) # porkchop_map[i,j]
    for i in range(n):
        for j in range(m):
            porkchop_map[i,j]=abs(data["vinf_arr"][i][j])+abs(data["vinf_dep"][i][j])
    
    """
    # Porchop Plot
    #fig, ax = plt.subplots()
    # fig = plt.figure()
    # plt.contour(dep_date, tof, porkchop_map, np.linspace(0,30,31), cmap="gnuplot") # arbitrary max contour level at dv = 30 km/s
    # fig.suptitle('Porkchop Plot for '+ data["object"]["fullname"])
    # plt.xlabel('$Date_{dep}$ (-)')
    # plt.ylabel('$ToF$ (d)')       
    # plt.colorbar()
    
    # locator = mdates.MonthLocator()
    # plt.gca().xaxis.set_major_locator(locator)
    
    # plt.gcf().autofmt_xdate()

    # date as xtick 
    # Major ticks every 6 months.
    # fmt_half_year = mdates.MonthLocator(interval=6)
    # ax.xaxis.set_major_locator(fmt_half_year)
    # Minor ticks every month.
    # fmt_month = mdates.MonthLocator()
    # ax.xaxis.set_minor_locator(fmt_month)
    # Text in the x axis will be displayed in 'YYYY-mm' format.
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    # Rotates and right aligns the x labels, and moves the bottom of the axes up to make room for them.
    #fig.autofmt_xdate()
    """

    return porkchop_map, dep_date, tof#, fig

##### DO NOT RUN NOW, TO CHECK THE USE
"""
# Mode T - Request list of small bodies that come closest (or within a prescribed distance) to a user-specified heliocentric 
# orbit (assumes two-body dynamics). Proxy for easiest-to-reach targets for an extended mission phase.
# example url
# https://ssd-api.jpl.nasa.gov/mdesign.api?ec=0.2056408220896557&qr=0.3074958016246215&tp=2459067.6508400026&
# om=48.30597718083336&w=29.18348714438387&in=7.003733902930839&jd0=2458849.5&jdf=2459132.5&maxout=100&maxdist=0.0010
"""
# MODE T
def get_close_approach_to_asteroid(orb_params,jd0,jdf,n_object_requested,distance_within):

    """
    # orb_params:            array containing the orbital parameters required to run the query
    #     ec:                eccentricity [>0]
    #     qr:                perihelion distance [>0]
    #     tp:                time of perihelion passage (JD)
    #     OM:                longitude of the ascending node (deg) [0,360]
    #     om:                argument of periapsis (deg) [0,360]
    #     incl:              inclination (deg) [0,180]
    # jd0:                   beginning of the requested time span (JD)
    # jdf:                   end of the requested time span (JD). Time span must not be longer than one year
    # n_object_requested:    maximum number of records to be returned
    # distance_within:       ignore objects with distance of closest approach greater than "distance_within" [>0, optional]
    """

    # Extraction of inputs
    ec = orb_params[1]
    qr = orb_params[2]
    tp = orb_params[3]
    OM = orb_params[4]
    om = orb_params[5]
    incl = orb_params[6]
    
    # Construction of the HTTP request
    url_base = 'https://ssd-api.jpl.nasa.gov/mdesign.api'
    url = f'{url_base}?ec={str(ec)}&qr={str(qr)}&tp={str(tp)}&om={str(OM)}&w={str(om)}&in={str(incl)}&jd0={str(jd0)}&jdf={str(jdf)}&maxout={str(n_object_requested)}&maxdist={str(distance_within)}'
    r = requests.get(url)
    data = r.json()

    return data

# MODE A
def get_accessible_sb(records_lim,optim_crit,years,sb_class,rdzvs,profile,ball_flag,lt_flag):

    """
    # If ballistic missions are requested, the API expects crit to be defined. 
    # If low-thrust missions are requested, the API expects profile and rdzvs to be defined.
    
    ## INPUT
    # lim:         number of records to be retrieved form the SBDB
    # crit:        optimality criterion for selecting ballistic missions: 
    #              1) min. departure V-infinity, 
    #              2) min. arrival V-infinity, 
    #              3) min. total delta-V, 
    #              4) min. tof + min. departure V-infinity, 
    #              5) min. tof + min. arrival V-infinity, 
    #              6) min. tof + min. total delta-V
    # year:        launch year or list of launch years for which optimal missions are to be retrieved 
    #              from the SBDB [current year + [0, 20]]
    # rdzvs        when requesting low-thrust missions, if rdzvs is true, only rendezvous missions are retrieved from the SBDB. 
    #              If false, flyby missions will be retrieved [boolean]
    # profile      when requesting low-thrust missions, profile maps to the spacecraft configuration: 
    #              1) Mid-size spacecraft, 2) smallsat
    
    ## OUTPUT BALLISTIC
    # name - small-body full name.
    # date0 - departure date (cal).
    # MJD0 - departure date (MJD).
    # datef - arrival date (cal).
    # MJDF - arrival date (MJD).
    # c3_dep - departure characteristic energy C3 (km^2/s^2).
    # vinf_dep - departure V-infinity (km/s).
    # vinf_arr - arrival V-infinity (km/s).
    # dv_tot - total delta-V (km/s).
    # tof - time of flight (d).
    # class - three-letter orbit class code.
    # H - absolute magnitude.
    # condition_code - orbit condition code.
    # neo - Near-Earth Object flag.
    # pha - Potentially-Hazardous Asteroid flag.
    # bin - binary flag.
    # pdes - designation.
    
    ## OUTPUT RENDEZ-VOUS
    # name - small-body full name.
    # date0 - departure date (cal).
    # MJD0 - departure date (MJD).
    # datef - arrival date (cal).
    # MJDF - arrival date (MJD).
    # c3_dep - departure characteristic energy C3 (km^2/s^2).
    # vinf_dep - departure V-infinity (km/s).
    # vinf_arr - arrival V-infinity (km/s).
    # mass_dep - departure mass.
    # mass_arr - arrival mass.
    # tof - time of flight (d).
    # nga - number of Earth gravity assists.
    # rdzvs - rendezvous flag.
    # class - three-letter orbit class code.
    # H - absolute magnitude.
    # condition_code - orbit condition code.
    # neo - Near-Earth Object flag.
    # pha - Potentially-Hazardous Asteroid flag.
    # bin - binary flag.
    # pdes - designation.
    """

    url_base = 'https://ssd-api.jpl.nasa.gov/mdesign.api'
    
    if (ball_flag==1 and lt_flag == 0): # ballistic profile mission requested
        # https://ssd-api.jpl.nasa.gov/mdesign.api?lim=200&crit=1&year=2025,2026,2027,2028,2029&sb_group=neo
        # Construction of the HTTP request
        url = f'{url_base}?lim={str(records_lim)}&crit={str(optim_crit)}&year={str(years)}&sb_group={str(sb_class)}'
    elif (ball_flag == 0 and lt_flag == 1): # low thrust profile mission requested
        # https://ssd-api.jpl.nasa.gov/mdesign.api?lim=200&rdzvs=true&profile=1&year=2025,2026,2027,2028,2029&sb_class=TJN
        # Construction of the HTTP request
        url = f'{url_base}?lim={str(records_lim)}&rdzvs={str(rdzvs)}&profile={str(profile)}&year={str(years)}&sb_class={str(sb_class)}'
        
    r = requests.get(url)
    data = r.json()
    
    # Elaboration of data
    accessible_valid_sb={}
    # Definition of limit condition to consider valid small body among the accessible ones
    mjd31Dec2028 = 62136
    mjd01Jan2048 = 69076
    lim_magnitude = 28
    lim_OCC = 7
    lim_c3_dep = 15 # km^2/s^2
    
    if (ball_flag==1 and lt_flag == 0): # ballistic profile mission requested
        j = 0
        for accessible_sb in range(int(data['md_constraints']['lim'])):
            if (int(data['data'][accessible_sb][2]) < mjd31Dec2028 and int(data['data'][accessible_sb][4]) < mjd01Jan2048 and float(data['data'][accessible_sb][11]) < lim_magnitude and int(data['data'][accessible_sb][12]) < lim_OCC and float(data['data'][accessible_sb][5]) < lim_c3_dep):
                valid_sb={"fullname": data['data'][accessible_sb][0],
                          "mjd0": data['data'][accessible_sb][2],
                          "mjdf": data['data'][accessible_sb][4],
                          "c3_dep": data['data'][accessible_sb][5],
                          "vinf_dep": data['data'][accessible_sb][6],
                          "vinf_arr": data['data'][accessible_sb][7],
                          "dv_tot": data['data'][accessible_sb][8],
                          "tof": data['data'][accessible_sb][9],
                          "class": data['data'][accessible_sb][10],
                          "H": data['data'][accessible_sb][11],
                          "condition_code": data['data'][accessible_sb][12],
                          }
                accessible_valid_sb[str(j)]=valid_sb
                j = j+1
    elif (ball_flag == 0 and lt_flag == 1): # low thrust profile mission requested
        j = 0
        for accessible_sb in range(int(data['md_constraints']['lim'])):
            if (int(data['data'][accessible_sb][2]) < mjd31Dec2028 and int(data['data'][accessible_sb][4]) < mjd01Jan2048 and float(data['data'][accessible_sb][14]) < lim_magnitude and int(data['data'][accessible_sb][15]) < lim_OCC and float(data['data'][accessible_sb][5]) < lim_c3_dep):
                valid_sb={"fullname": data['data'][accessible_sb][0],
                          "mjd0": data['data'][accessible_sb][2],
                          "mjdf": data['data'][accessible_sb][4],
                          "c3_dep": data['data'][accessible_sb][5],
                          "vinf_dep": data['data'][accessible_sb][6],
                          "vinf_arr": data['data'][accessible_sb][7],
                          "mass_dep": data['data'][accessible_sb][8],
                          "mass_arr": data['data'][accessible_sb][9],
                          "tof": data['data'][accessible_sb][10],
                          "nga": data['data'][accessible_sb][11],
                          "rdzvs": data['data'][accessible_sb][12],
                          "class": data['data'][accessible_sb][13],
                          "H": data['data'][accessible_sb][14],
                          "condition_code": data['data'][accessible_sb][15],
                          }
                accessible_valid_sb[str(j)]=valid_sb;
                j = j+1        

    return accessible_valid_sb

def get_min_dv_accessible_sb(accessible_valid_sb):

    # Find the best Mission Profile among the valid accessible small bodies
    dv = np.zeros(len(accessible_valid_sb.keys()))
    for profile in range(len(accessible_valid_sb)):
        dv[profile] = accessible_valid_sb[str(profile)]['dv_tot']

    index = np.linspace(0,len(dv)-1,len(dv))
    mask = [dv == np.min(dv)]
    idx_min = index[tuple(mask)]
    if type(idx_min) == float:
        accessible_sb_min_dv = accessible_valid_sb[str(int(idx_min))]
    else:
        #val, idx_min = min((val, idx) for (idx, val) in enumerate(dv))
        accessible_sb_min_dv = accessible_valid_sb[str(int(idx_min[0]))] # NOTE: there could be more than one best solution, here i took arbitrarly the first
    
    # Plot of the mission profiles and highlight the best one
    fig = plt.figure()
    plt.plot(dv, "*")
    fig.suptitle('Accessible Small Bodies Mission dv, and the minimum is ' + accessible_sb_min_dv["fullname"])
    plt.xlabel('$idx$ (-)')
    plt.ylabel('$dv$ (km/s)')
    plt.plot(int(idx_min[0]), dv[int(idx_min[0])], "r+")
    
    return accessible_sb_min_dv, fig

# JPL SBDB 
def get_dict(name_list):

    """
    # Return information for each bodies in name_list from JPL Small Body Database in form of dictionair
    # INPUT
    # name_list      list [name1, name2, ..., nameN]
    # OUTPUT
    # dict_bodies    dict with the following structure
    # str(name): "fullname"
    #         "our_id"                is the id inside the new dict, from 0 to N-1 where N is the name_list length
    #         "neo_flag"
    #         "orbit_class"
    #         "pha_flag"
    #         "object_kind"
    #         "moid"
    #         "orbita_elements"
    #         "condition_code"
    #         "rms"
    #         "orbit_comment"
    #         "magn_radius_flag"      is "H" if the next parameter is the magnitude, is "D" if the next parameter is the diameter
    #         "H" or "D"
    #         "spectral_category_flag"     T, S or 0
    #         "spectral_category"
    #         "N_obs"                      Number of observations
    #         "obs_span"                   Time between first and last observation
    #         "impacts": str(impact id):   'width'  
    #                                      'energy'
    #                                      'stretch'
    #                                      'ip'
    #                                      'dt'
    #                                      'date'
    #                                      'sigma_lov'
    #                                      'h'
    #                                      'mass'
    #                                      'v_inf'
    #                                      'sigma_imp'
    #                                      'method'
    #                                      'ts'
    #                                      'diam'
    #                                      'dist'
    #                                      'v_imp'
    #                                      'ps'
    """

    our_id=0
    dict_bodies={}

    for name in name_list:
        sbdb = SBDB.query(name, neo_only=True, full_precision=True, phys=True, virtual_impactor=True)
        if sbdb["object"]["kind"]!='cn' or sbdb["object"]["kind"]!='cu' :
            asteroid={"fullname": sbdb["object"]["fullname"],# TODO vogliamo fare un check?
              "our_id":our_id,
              "neo_flag": sbdb["object"]["neo"],
              "orbit_class":sbdb["object"]["orbit_class"]["code"],
              "pha_flag":sbdb["object"]["pha"],
              "object_kind":sbdb["object"]["kind"], #an asteroid numbered au unbered asteroid (cn, cu for comet)
              "moid": sbdb["orbit"]["moid"],
              "orbital_elements":sbdb["orbit"]['elements'],
              "condition_code": sbdb["orbit"]["condition_code"], #OCC
              "rms": sbdb["orbit"]["rms"],
              "orbit_comment":sbdb["orbit"]["comment"],
             }
            try:
                asteroid["magn_radius_flag"]='H'
                asteroid["H"]=sbdb['phys_par']['H']
            except:
                asteroid["magn_radius_flag"]='D'
                asteroid["D"]=sbdb['phys_par']['diameter']
            asteroid["N_obs"]=sbdb['orbit']['n_obs_used']
            asteroid["obs_span"]=sbdb['orbit']['data_arc']
            
            asteroid["impacts"]={}
            flag_bool=1
            if 'phys_par' in sbdb.keys():
                spect_flag=0
                if 'spec_T' in sbdb['phys_par'].keys():
                    asteroid["spectral_category_flag"]='T'
                    asteroid["spectral_category"]=sbdb['phys_par']['spec_T']
                    spect_flag=1
                if 'spec_B' in sbdb['phys_par'].keys():
                    asteroid["spectral_category_flag"]='B'
                    asteroid["spectral_category"]=sbdb['phys_par']['spec_B']
                    spect_flag=1
                if spect_flag==0:
                    asteroid["spectral_category_flag"]='0'
            else:
                asteroid["spectral_category"]='0'
                
            if 'ip' in sbdb["vi_data"]:
                n_imp=len(sbdb["vi_data"]['ip'])
                for key in sbdb["vi_data"].keys():
                    if flag_bool==1:
                        for i in range(0,n_imp):
                            asteroid["impacts"][str(i)]={}
                        flag_bool=0
                    for i in range(0,n_imp): 
                        try:
                            if isinstance(sbdb["vi_data"][key],str):
                                asteroid["impacts"][str(i)][key]=sbdb["vi_data"][key];
                            else:
                                asteroid["impacts"][str(i)][key]=sbdb["vi_data"][key][i];
                        except:
                            print(name+" could raise error in importing virtual impact data") #this exception is raised if only one impact is present
            dict_bodies[name]=asteroid
            our_id=our_id+1
            flag_bool=1
            del asteroid

    return dict_bodies
        
def extract_esa_name_from_file(file_name):

    f = open(file_name, "r")
    #f = open('NEO_API_py/' + file_name, "r")
    line=f.readline()
    line=f.readline()
    line=f.readline()
    line=f.readline()
    counter=0
    esa_names=[]
    for line in f:
        word=""
        for c in line:
            if c==' ':
                break
            else:
                word=word+c
        esa_names.append(word);  

    return esa_names

def get_sentry_risk_list():

    url = 'https://ssd-api.jpl.nasa.gov/sentry.api'
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/81.0.4044.138 Safari/537.36'}

    r = requests.get(url, headers=headers)
    data = r.json()
    sentry_risk_names=[]
    for i in range(0,len(data['data'])):
        name_=""
        name=data['data'][i]['des']
        for c in name:
            if c==' ':
                a=0
            else:
                name_=name_+c
        sentry_risk_names.append(name_)

    return sentry_risk_names

def merge_risk_lists(esa,sentry):
    risk_list=sentry
    for risk_name in esa:
        if risk_name not in sentry:
            risk_list.append(risk_name)
    return risk_list

def refined_selection(dict_risk_list):

    # This function return a dictionair of asteroid satisfyng the requirements

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
            max_P=-100
            date_flag=0
            for imp_id in dict_risk_list[key]['impacts'].keys():
                word=''
                for c in dict_risk_list[key]['impacts'][imp_id]['date']:
                    #pprint(c)
                    if c=='-':
                        break
                    else:
                        word=word+c
                if int(word)<2048 and int(word)>2026:
                    date_flag=1
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
        print(refined_selected[ind])
        print('Max Palermo Scale:' + str(PS_list[ind]))
        print('OCC:' + str(dict_risk_list[refined_selected[ind]]['condition_code']))
        refined_dict[refined_selected[ind]]={}
        refined_dict[refined_selected[ind]]=dict_risk_list[refined_selected[ind]]

    return(refined_dict, refined_selected)

# DB EXPLORATION TOOL
def MOID_H(dict_risk_list):

    #Earth minimum orbit instersection distance and magnitude plotting
    H_=[]
    MOID_=[]
    for key in dict_risk_list:
        try:
            H_.append(float(dict_risk_list[key]['H']))
            MOID_.append(float(dict_risk_list[key]['moid'].scale))
        except:
            print(key+' does not have magnitude info')
    x = MOID_
    y = H_

    """
    # fig, ax = plt.subplots()
    # ax.plot(x,y, marker='o', linewidth=0)
    # start, end = ax.get_ylim()
    # ax.yaxis.set_ticks(np.arange(-14, 306, 20))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    # plt.show()
    """

    return x,y

def H_OCC(dict_risk_list):

    #Magnitude and Orbital Condition Code plotting
    H_=[]
    OCC_=[]
    for key in dict_risk_list:
        try:
            H_.append(float(dict_risk_list[key]['H']))
            OCC_.append(int(dict_risk_list[key]['condition_code']))
        except:
            print(key+' does not have magnitude info')
    x = H_
    y = OCC_

    """
    # fig, ax = plt.subplots()
    # ax.plot(x,y, marker='o', linewidth=0)
    # start, end = ax.get_ylim()
    # ax.yaxis.set_ticks(np.arange(-14, 306, 20))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    # plt.show()
    """
    return x,y

def get_df_for_sbdb_visualization(dict_risk_list):
    moid=[]; occ=[]; H=[]; worse_impact_ps=[];
    idx = 0
    arbitrary_albedo = 0.14 #between 0.05 and 0.25, but most of the asteroids are on higher ranges
    for el in dict_risk_list:
        moid.append(float(dict_risk_list[str(el)]['moid'].scale))
        occ.append(int(dict_risk_list[str(el)]['condition_code']))
        try:
            H.append(float(dict_risk_list[str(el)]['H']))
        except:
            print(el +' does not have magnitude info')
            #http://www.physics.sfasu.edu/astro/asteroids/sizemagnitude.html
            H.append(float(-np.log10(dict_risk_list[str(el)]['D'].scale*np.sqrt(arbitrary_albedo)/1329)*5));

        worse_impact_ps_lim = -20
        worse_impact_ps.append(float(worse_impact_ps_lim))

        for i in range(len(dict_risk_list[str(el)]['impacts'])):
            if (dict_risk_list[str(el)]['impacts'][str(i)]['ps'] > worse_impact_ps_lim):
                worse_impact_ps.pop(idx)
                worse_impact_ps.insert(idx, float(dict_risk_list[str(el)]['impacts'][str(i)]['ps']))
                
                worse_impact_ps_lim = float(dict_risk_list[str(el)]['impacts'][str(i)]['ps'])     
                
        idx = idx + 1

    kinda_dict_physical_properties = {'moid': array("f",moid), 'occ': array("i",occ), 
                                    'H':  array("f",H),'worse_impact_ps': array("f",worse_impact_ps)}
    df_physical_properties = pd.DataFrame(data=kinda_dict_physical_properties)

    physical_properties_name = kinda_dict_physical_properties.keys()

    return df_physical_properties, physical_properties_name

def bi_impulsive_mission(refined_selected, mjd0, duration, min_tof, max_tof, step_size):
    
    """
    mjd0 MJD2000 of departure date
    duration
    min_tof minimum time of flight
    max_tof maximum time of flight
    step_size step size
    """

    refined_selected_MD={}
    for name in refined_selected:
        refined_selected_MD[name]={}
        refined_selected_MD[name]['missions'], refined_selected_MD[name]['porkchop_dv'], refined_selected_MD[name]['dep_date'], refined_selected_MD[name]['tof'],  refined_selected_MD[name]['mp_min_dv'] = \
            get_mission_profiles(name,mjd0,duration,min_tof,max_tof,step_size) #removed refined_selected_MD[name]['pc_plot'], refined_selected_MD[name]['mp_dv_plot'] check order
    
    return refined_selected_MD

# NASA JPL Horizons
def get_horizons_ephemerides(name,pov,epoch_start,epoch_stop,step_size,type_elements):
    
    # step: step size, [10m, 1d, 1y]
    
    if pov.lower() == 'sun':
        loc = '500@10' # position relative to the sun
    elif pov.lower() == 'goldstone':
        loc = '257' # from goldstone
    elif pov.lower() == 'maunakea':
        loc = '568' # maunakea
    else:
        print('Not Valid Location Point Of View')
    
    # Process to get homogeneity from main script full name '2012QD8' to a valid name for Horizon call '2012 QD8'
    if len(re.findall('([0-9])', name)) <= 4: # 4 is the min numbers in every name, the date year of discovery
        r = re.compile("([0-9]+)([a-zA-Z]+)").match(name)
        k1 = r.group(1) # the date of the name
        k2 = r.group(2) # the code of the date
        valid_name = k1 + " " + k2 
    else:
        r = re.compile("([0-9]+)([a-zA-Z]+)([0-9]+)").match(name)
        k1 = r.group(1) # the date of the name
        k2 = r.group(2) # the code of the date
        k3 = r.group(3) # id after the letters
        valid_name = k1 + " " + k2 + k3
    
    obj = Horizons(id=valid_name, 
               location=loc, 
               epochs={'start': epoch_start, 'stop':epoch_stop,
                       'step': step_size})
    
    if type_elements.lower() == 'vectors':
        data_output = obj.vectors() # vectorial elements
    elif type_elements.lower() == 'ephemerides':
        data_output = obj.ephemerides()

    len_rows = len(data_output)
    len_cols = 6 # 3 positions 'x','y','z', and 3 velocities 'vx', 'vy', 'vz'
    idx_x = 5 # 'x' is at position 5 in the table (starting from 0)
    data =  np.zeros([len_rows,len_cols]) 
    for row in range(len_rows):
        for col in range(6): 
            idx_col_in_table = idx_x + col # because the 'x' start at 6th col, going up till the 12th that is 'vz'
            data[row,col] = data_output[row][idx_col_in_table]      
            
    return data

def get_earth_ephemerides(epoch_start,epoch_stop,step_size,type_elements):
    
    # step: step size, [10m, 1d, 1y]

    obj = Horizons(id = 'Geocenter', 
               location = '500@10', 
               epochs = {'start': epoch_start, 'stop':epoch_stop,
                       'step': step_size},
               id_type = 'majorbody')
    
    if type_elements.lower() == 'vectors':
        data_output = obj.vectors() # vectorial elements
    elif type_elements.lower() == 'ephemerides':
        data_output = obj.ephemerides()
    
    len_rows = len(data_output)
    len_cols = 6 # 3 positions 'x','y','z', and 3 velocities 'vx', 'vy', 'vz'
    idx_x = 3 # 'x' is at position 3 in the table
    data =  np.zeros([len_rows,len_cols]) 
    for row in range(len_rows):
        for col in range(6): 
            idx_col_in_table = idx_x + col # because the 'x' start at 3rd col, going up till the 9th that is 'vz'
            data[row,col] = data_output[row][idx_col_in_table]

    return data

