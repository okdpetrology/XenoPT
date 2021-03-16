#This project was made in December 2020 by Khalil Droubi for UW-Madison, Fall2 2020, CS319 class AND for research purposes. 

#Unless given explicit permission by Khalil, you do NOT have permission to share this code with anybody.

#imports
import math
import pandas as pd
import os
import csv
#import copy
#import sqlite3
import warnings

# Define Molar mass values in dictionary
Oxide_Weight_conversions = {
    "Si": 60.0843, 
    "Ti": 79.8988, 
    "Al": 101.96128, 
    "Cr": 151.99, 
    "Fe": 71.8464, 
    "Mn":70.9374, 
    "Ni": 74.6994, 
    "Mg": 40.3044, 
    "Ca": 56.0794, 
    "Na": 61.97894,
    "K": 94.196}

owc = Oxide_Weight_conversions

def moles_oxygen(dataframe):

    xeno_master_work = {}
    for idx in range(len(dataframe)): 
        xeno_work = {}
        for row in dataframe.columns:

            for element in Oxide_Weight_conversions:
                if element in row:
                    #print(element + ": " + row + ": " + str(xeno.loc[idx, row]))
                    rename = row.split("_")[0] + "_" + element
                    if element == "Si":
                        #print(rename)
                        xeno_work[rename] = dataframe.loc[idx, row] / owc["Si"] * 2
                    if element == "Ti":
                        xeno_work[rename] = dataframe.loc[idx, row] / owc["Ti"] * 2    
                    if element == "Al":
                        xeno_work[rename] = dataframe.loc[idx, row] / owc["Al"] * 3 
                    if element == "Cr":
                        xeno_work[rename] = dataframe.loc[idx, row] / owc["Cr"] * 3 
                    if element == "Fe":
                        xeno_work[rename] = dataframe.loc[idx, row] / owc["Fe"] 
                    if element == "Mn":
                        xeno_work[rename] = dataframe.loc[idx, row] / owc["Mn"] 
                    if element == "Ni":
                        xeno_work[rename] = dataframe.loc[idx, row] / owc["Ni"] 
                    if element == "Mg":
                        xeno_work[rename] = dataframe.loc[idx, row] / owc["Mg"] 
                    if element == "Ca":
                        xeno_work[rename] = dataframe.loc[idx, row] / owc["Ca"]
                    if element == "Na":
                        xeno_work[rename] = dataframe.loc[idx, row] / owc["Na"] 
                    if element == "K":
                        xeno_work[rename] = dataframe.loc[idx, row] / owc["K"]

        name = dataframe.loc[idx, "SampleID"]   
        xeno_master_work[name] = xeno_work
        
    return xeno_master_work

def sum_and_factor(dd): 
    #Input the dict of dicts (dd) from moles_oxygen definition
    mineral_sum = {}
    mineral_factor = {}
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="invalid value encountered in double_scalars")
        warnings.filterwarnings("ignore", message="divide by zero encountered in double_scalars")
        for sample in dd:
            sample_sum = {"cpx": 0, "opx": 0, "olv": 0, "gt": 0}
            sample_factor = {}
            for entry in dd[sample]:
                #print(entry)
                #print(df[sample][entry])
                if "cpx" in entry:
                    sample_sum["cpx"] = sample_sum["cpx"] + dd[sample][entry]
                if "opx" in entry:
                    sample_sum["opx"] = sample_sum["opx"] + dd[sample][entry]
                if "olv" in entry:
                    sample_sum["olv"] = sample_sum["olv"] + dd[sample][entry]
                if "gt" in entry:
                    sample_sum["gt"] = sample_sum["gt"] + dd[sample][entry]

            sample_factor["cpx"] = 1 / sample_sum["cpx"] * 6
            sample_factor["opx"] = 1 / sample_sum["opx"] * 6
            sample_factor["olv"] = 1 / sample_sum["olv"] * 4
            sample_factor["gt"]  = 1 / sample_sum["gt"]  * 12

            mineral_sum[sample] = sample_sum
            mineral_factor[sample] = sample_factor
    return mineral_sum, mineral_factor

def formula_units(dd, sample_factor):
    #Input the dict of dicts (dd) from moles_oxygen definition
    sample_formula_units = {}
    
    for sample in dd:
        mineral_formula_units = {}
        
        for entry in dd[sample]:
                   #print(entry)
            if "cpx" in entry:
                working_factor = sample_factor[sample]["cpx"]
            if "opx" in entry:
                working_factor = sample_factor[sample]["opx"]
            if "olv" in entry:
                working_factor = sample_factor[sample]["olv"]
            if "gt" in entry:
                working_factor = sample_factor[sample]["gt"]

            wf = working_factor
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="invalid value encountered in double_scalars")
                warnings.filterwarnings("ignore", message="divide by zero encountered in double_scalars")
                for element in owc:
                    if element in entry:

                        if element == "Si":
                            mineral_formula_units[entry] = dd[sample][entry] * wf / 2

                        if element == "Ti":
                            mineral_formula_units[entry] = dd[sample][entry] * wf / 2

                        if element == "Al":
                            mineral_formula_units[entry] = dd[sample][entry] * wf * 2/3

                        if element == "Cr":
                            mineral_formula_units[entry] = dd[sample][entry] * wf * 2/3

                        if element == "Fe":
                            mineral_formula_units[entry] = dd[sample][entry] * wf

                        if element == "Mn":
                            mineral_formula_units[entry] = dd[sample][entry] * wf

                        if element == "Ni":
                            mineral_formula_units[entry] = dd[sample][entry] * wf

                        if element == "Mg":
                            mineral_formula_units[entry] = dd[sample][entry] * wf

                        if element == "Ca":
                            mineral_formula_units[entry] = dd[sample][entry] * wf

                        if element == "Na":
                            mineral_formula_units[entry] = dd[sample][entry] * wf * 2
                        if element == "K":
                            mineral_formula_units[entry] = dd[sample][entry] * wf * 2
                           
        sample_formula_units[sample] = mineral_formula_units
        
    return sample_formula_units

def formula_units_fast(dataframe):
    dd = moles_oxygen(dataframe)
    #dd = dict of dicts
    x = list(sum_and_factor(dd))
    sample_factor = x[1]
    
    dc = formula_units(dd, sample_factor)
    
    return dc

#Site Occupancy


def site_occupancy(wfu, ferric_choice = 'stoichiometry'):
    '''
    Possible inputs for ferric_choice are:
    (default input) stoichiometry = Calculate Fe3+ in cpx and opx based on cation stoichiometry.
    kd_assumption =  Estimates of maximum Fe3+ as ~25% total iron and ~15% total iron, for cpx and opx, according to the range of data in Dyar et al. (1989) and Luth & Canil (1992).
    [cpx_%Fe3+, opx_%Fe3+] = input a list containing values from 0 to 1 to manually set what percentage of total iron is Fe3+ in each phase.
    '''
    site_occupancy = {}
    fc = ferric_choice
    for sample in wfu:
        min_site = {}
        ms = min_site

        #clinopyroxene and orthopyroxene site occupancy

        ms["cpx_Mg/(Mg+Fe)"] = wfu[sample]["cpx_Mg"] / (wfu[sample]["cpx_Mg"] + wfu[sample]["cpx_Fe"])
        ms["opx_Mg/(Mg+Fe)"] = wfu[sample]["opx_Mg"] / (wfu[sample]["opx_Mg"] + wfu[sample]["opx_Fe"])

        ms["cpx_Mg+Fe_M1"] = 3 - wfu[sample]["cpx_Al"] -  wfu[sample]["cpx_Si"] -  wfu[sample]["cpx_Cr"] - wfu[sample]["cpx_Ti"]
        ms["opx_Mg+Fe_M1"] = 3 - wfu[sample]["opx_Al"] -  wfu[sample]["opx_Si"] -  wfu[sample]["opx_Cr"] - wfu[sample]["opx_Ti"]
        ms["cpx_Mg+Fe_M2"] = 1 - wfu[sample]["cpx_Mn"] -  wfu[sample]["cpx_Ca"] -  wfu[sample]["cpx_Na"] - wfu[sample]["cpx_K"]
        ms["opx_Mg+Fe_M2"] = 1 - wfu[sample]["opx_Mn"] -  wfu[sample]["opx_Ca"] -  wfu[sample]["opx_Na"] - wfu[sample]["opx_K"]

        ms["cpx_Mg_M2"] = ms["cpx_Mg+Fe_M2"] * ms["cpx_Mg/(Mg+Fe)"]
        ms["opx_Mg_M2"] = ms["opx_Mg+Fe_M2"] * ms["opx_Mg/(Mg+Fe)"]

        ms["cpx_a_Mg"] = ms["cpx_Mg_M2"] * (wfu[sample]["cpx_Mg"] - ms["cpx_Mg_M2"])
        ms["opx_a_Mg"] = ms["opx_Mg_M2"] * (wfu[sample]["opx_Mg"] - ms["opx_Mg_M2"])

        ms["opx_XM1_Al_HG"] = wfu[sample]["opx_Al"] / 2
        ms["cpx_XM1_Al_NG"] = (wfu[sample]["cpx_Al"] - wfu[sample]["cpx_Cr"] - (2 * wfu[sample]["cpx_Ti"]) + wfu[sample]["cpx_Na"]) / 2
        ms["opx_XM1_Al_NG"] = (wfu[sample]["opx_Al"] - wfu[sample]["opx_Cr"] - (2 * wfu[sample]["opx_Ti"]) + wfu[sample]["opx_Na"]) / 2

        ms["opx_XM1_Fe,Mg"] = 1 - ms["opx_XM1_Al_NG"] - wfu[sample]["opx_Cr"] - wfu[sample]["opx_Ti"]
        ms["opx_XM1_Mg"] = ms["opx_XM1_Fe,Mg"] * ms["opx_Mg/(Mg+Fe)"]
        ms["opx_XM1_Fe"] = ms["opx_XM1_Fe,Mg"] * (wfu[sample]["opx_Fe"] / (wfu[sample]["opx_Fe"] + wfu[sample]["opx_Mg"]))
        ms["opx_XM1_Cr"] = 1 - ms["opx_XM1_Fe,Mg"] - ms["opx_XM1_Al_NG"] - wfu[sample]["opx_Ti"]

        ms["cpx_Fe#"] = wfu[sample]["cpx_Fe"] / (wfu[sample]["cpx_Fe"] + wfu[sample]["cpx_Mg"])
        ms["opx_Fe#"] = wfu[sample]["opx_Fe"] / (wfu[sample]["opx_Fe"] + wfu[sample]["opx_Mg"])

        ms["cpx_Fe/Mg"] = wfu[sample]["cpx_Fe"] / wfu[sample]["cpx_Mg"]
        ms["opx_Fe/Mg"] = wfu[sample]["opx_Fe"] / wfu[sample]["opx_Mg"]

        # Calculating Fe3+ based on cation stoichiometry
        # I don't agree with this for pyroxenes, but we are going to start out with this calculation, 
        # and figure out how to manually change it later

        # Working on Fe3+ fix, but not finished

        if fc == 'stoichiometry':
            ms["cpx_Fe3+"] = (2 - wfu[sample]["cpx_Si"]) + wfu[sample]["cpx_Na"] - (wfu[sample]["cpx_Al"] - 2 + wfu[sample]["cpx_Si"] + wfu[sample]["cpx_Cr"] + (2 * wfu[sample]["cpx_Ti"]))
            ms["opx_Fe3+"] = (2 - wfu[sample]["opx_Si"]) + wfu[sample]["opx_Na"] - (wfu[sample]["opx_Al"] - 2 + wfu[sample]["opx_Si"] + wfu[sample]["opx_Cr"] + (2 * wfu[sample]["opx_Ti"]))

        if fc == 'kd_assumption':
            ## Estimates of maximum Fe3+ as ~25% total iron and ~15% total iron, for clinopyroxene and orthopyroxene,
            # according to the range of data in Dyar et al. (1989) and Luth & Canil (1992).
            ms["cpx_Fe3+"] = wfu[sample]["cpx_Fe"] * 0.25
            ms["opx_Fe3+"] = wfu[sample]["opx_Fe"] * 0.15
        # Future if statements should allow for manual set of cpx or opx Fe3+ = X% of Fe(total)
        if type(fc) == "list":
            ms["cpx_Fe3+"] = wfu[sample]["cpx_Fe"] * fc[0]
            ms["opx_Fe3+"] = wfu[sample]["opx_Fe"] * fc[1]
        else:
            pass
            #ms["cpx_Fe3+"] = pass
            #ms["opx_Fe3+"] = pass

        ms["cpx_%Fe3+"] = 100 * ms["cpx_Fe3+"] / wfu[sample]["cpx_Fe"]
        ms["opx_%Fe3+"] = 100 * ms["opx_Fe3+"] / wfu[sample]["opx_Fe"]   

        # We want to be able to  change the %Fe3+ to between 0 and 24%, realistically.

        ms["cpx_Fe2+/Mg"] = (wfu[sample]["cpx_Fe"] - ms["cpx_Fe3+"]) / wfu[sample]["cpx_Mg"]
        ms["opx_Fe2+/Mg"] = (wfu[sample]["opx_Fe"] - ms["opx_Fe3+"]) / wfu[sample]["opx_Mg"]
        ms["cpx_Jadeite"] = wfu[sample]["cpx_Na"] - wfu[sample]["cpx_Cr"] - ms["cpx_Fe3+"] - (2 * wfu[sample]["cpx_Ti"])
        ms["opx_Jadeite"] = wfu[sample]["opx_Na"] - wfu[sample]["opx_Cr"] - ms["opx_Fe3+"] - (2 * wfu[sample]["opx_Ti"])

        ms["cpx_XAl_M1"] = (wfu[sample]["cpx_Al"] + abs(ms["cpx_Jadeite"])) / 2
        ms["opx_XAl_M1"] = (wfu[sample]["opx_Al"] + ms["opx_Jadeite"]) / 2

        ms["cpx_Xmf_M1"] = 1 - ((wfu[sample]["cpx_Al"] + ms["cpx_Jadeite"]) / 2) -  wfu[sample]["cpx_Cr"] - wfu[sample]["cpx_Ti"]
        ms["opx_Xmf_M1"] = 1 - ms["opx_XAl_M1"] - wfu[sample]["opx_Cr"] - ms["opx_Fe3+"] - wfu[sample]["opx_Ti"]
        ms["cpx_Xmf_M2"] = 1 - wfu[sample]["cpx_Ca"] - wfu[sample]["cpx_Na"] - wfu[sample]["cpx_Mn"]
        ms["opx_Xmf_M2"] = 1 - wfu[sample]["opx_Ca"] - wfu[sample]["opx_Na"] - wfu[sample]["opx_Mn"]

        ms["cpx_Xmg_M1"] = ms["cpx_Xmf_M1"] * ms["cpx_Mg/(Mg+Fe)"]
        ms["opx_Xmg_M1"] = ms["opx_Mg/(Mg+Fe)"] * (1 - ((wfu[sample]["opx_Al"] + ms["opx_Jadeite"]) / 2) - wfu[sample]["opx_Cr"] - wfu[sample]["opx_Ti"])
        ms["cpx_Xmg_M2"] = ms["cpx_Xmf_M2"] * ms["cpx_Mg/(Mg+Fe)"]
        ms["opx_Xmg_M2"] = ms["opx_Mg/(Mg+Fe)"] * (1 - wfu[sample]["opx_Ca"] - wfu[sample]["opx_Na"] - wfu[sample]["opx_Mn"])
        ms["cpx_Xfe_M1"] = (ms["cpx_Xmg_M1"] - (ms["cpx_Xmg_M1"] * ms["cpx_Mg/(Mg+Fe)"])) / ms["cpx_Mg/(Mg+Fe)"]
        ms["opx_Xfe_M1"] = (ms["opx_Xmg_M1"] - (ms["opx_Xmg_M1"] * ms["opx_Mg/(Mg+Fe)"])) / ms["opx_Mg/(Mg+Fe)"]

        ms["cpx_a(en)"] = (1 - wfu[sample]["cpx_Ca"] - wfu[sample]["cpx_Na"] - wfu[sample]["cpx_K"]) * ( 1 - (0.5 * (wfu[sample]["cpx_Al"] + wfu[sample]["cpx_Cr"] + wfu[sample]["cpx_Na"] + wfu[sample]["cpx_K"])))                                                                                               

        #garnet site occupancy

        ms["gt_XCa"] = wfu[sample]["gt_Ca"] / (wfu[sample]["gt_Ca"] + wfu[sample]["gt_Mg"] + wfu[sample]["gt_Fe"] + wfu[sample]["gt_Mn"])
        ms["gt_XAl"] = wfu[sample]["gt_Al"] / (wfu[sample]["gt_Al"] + wfu[sample]["gt_Cr"])
        ms["gt_XCr"] = wfu[sample]["gt_Cr"] / (wfu[sample]["gt_Al"] + wfu[sample]["gt_Cr"])
        ms["gt_Fe#"] = wfu[sample]["gt_Fe"] / (wfu[sample]["gt_Fe"] + wfu[sample]["gt_Mg"])
        ms["gt_XFe"] = wfu[sample]["gt_Fe"] / (wfu[sample]["gt_Fe"] + wfu[sample]["gt_Mn"] + wfu[sample]["gt_Mg"] + wfu[sample]["gt_Ca"])
        ms["gt_XMg"] = wfu[sample]["gt_Mg"] / (wfu[sample]["gt_Fe"] + wfu[sample]["gt_Mn"] + wfu[sample]["gt_Mg"] + wfu[sample]["gt_Ca"])
        ms["gt_XMn"] = wfu[sample]["gt_Mn"] / (wfu[sample]["gt_Fe"] + wfu[sample]["gt_Mn"] + wfu[sample]["gt_Mg"] + wfu[sample]["gt_Ca"])
        ms["gt_Fe/Mg"] = wfu[sample]["gt_Fe"] / wfu[sample]["gt_Mg"]

        #olivine site occupancy

        ms["olv_Fe#"] = wfu[sample]["olv_Fe"] / (wfu[sample]["olv_Fe"] + wfu[sample]["olv_Mg"])
        ms["olv_Fe/Mg"] = wfu[sample]["olv_Fe"] / wfu[sample]["olv_Mg"]


        site_occupancy[sample] = min_site

    return site_occupancy

# Beyer site occupancy

# note that Fe3+ will affect these results, so need to figure out how to manually change Fe3+ in wfu

def beyer_paramaters(wfu, so):
    '''
    Calculated parameters from the Beyer 2016 ResearchGate spreadsheet:

    https://www.researchgate.net/publication/306038208_Update_2_Beyer_et_al_2015_grt-cpx_barometer_for_eclogites_Iterative_calculation_of_P_and_T
    '''
    site_occupancy = {}

    for sample in wfu:
        min_site = {}
        ms = min_site
        
        #garnet site occupancy
        
        ms['gt_normFe'] = 3 * wfu[sample]['gt_Fe'] / (wfu[sample]['gt_Fe'] + wfu[sample]['gt_Mg'] + wfu[sample]['gt_Ca'])
        ms['gt_normMg'] = 3 * wfu[sample]['gt_Mg'] / (wfu[sample]['gt_Fe'] + wfu[sample]['gt_Mg'] + wfu[sample]['gt_Ca'])
        ms['gt_normCa'] = 3 * wfu[sample]['gt_Ca'] / (wfu[sample]['gt_Fe'] + wfu[sample]['gt_Mg'] + wfu[sample]['gt_Ca'])
        ms['gt_catTotal'] = ms['gt_normFe'] + ms['gt_normMg'] + ms['gt_normCa']
        
        ms['gt_XFe'] = wfu[sample]['gt_Fe'] / (wfu[sample]['gt_Fe'] + wfu[sample]['gt_Mg'] + wfu[sample]['gt_Ca'])
        ms['gt_XMg'] = wfu[sample]['gt_Mg'] / (wfu[sample]['gt_Fe'] + wfu[sample]['gt_Mg'] + wfu[sample]['gt_Ca'])
        ms['gt_XCa'] = wfu[sample]['gt_Ca'] / (wfu[sample]['gt_Fe'] + wfu[sample]['gt_Mg'] + wfu[sample]['gt_Ca'])

    
        #clinopyroxene site occupancy
        
        ms['cpx_Al_IV'] = 2 - wfu[sample]['cpx_Si']
        ms['cpx_Fe+Mg_M1'] = 1 - (wfu[sample]['cpx_Al'] - ms['cpx_Al_IV']) - so[sample]['cpx_Fe3+'] - wfu[sample]['cpx_Ti'] 
        ms['cpx_XFe2+_M1'] = wfu[sample]['cpx_Fe'] / (wfu[sample]['cpx_Fe'] + wfu[sample]['cpx_Mg']) * ms['cpx_Fe+Mg_M1']
        ms['cpx_XFe3+_M1'] = so[sample]['cpx_Fe3+']
        ms['cpx_XMg_M1'] = wfu[sample]['cpx_Mg'] / (wfu[sample]['cpx_Fe'] + wfu[sample]['cpx_Mg']) * ms['cpx_Fe+Mg_M1']
        ms['cpx_XAl_M1'] = wfu[sample]['cpx_Al'] - ms['cpx_Al_IV']
        
        ms['cpx_Fe+Mg_M2'] = wfu[sample]['cpx_Fe'] + wfu[sample]['cpx_Mg'] - ms['cpx_Fe+Mg_M1']
        ms['cpx_XCa_M2'] = wfu[sample]['cpx_Ca']
        ms['cpx_XNa_M2'] = wfu[sample]['cpx_Na']
        ms['cpx_XFe_M2'] = wfu[sample]['cpx_Fe'] / (wfu[sample]['cpx_Fe'] + wfu[sample]['cpx_Mg']) * ms['cpx_Fe+Mg_M2']
        ms['cpx_XMg_M2'] = wfu[sample]['cpx_Mg'] / (wfu[sample]['cpx_Fe'] + wfu[sample]['cpx_Mg']) * ms['cpx_Fe+Mg_M2']
        ms['cpx_XAl_Esk_M2'] = (wfu[sample]['cpx_Al'] - wfu[sample]['cpx_Na']) - (2 *(2 - wfu[sample]['cpx_Si']))
            #XAl_Esk is weird...
        ms['cpx_Ca+Na_M2'] = ms['cpx_XCa_M2'] + ms['cpx_XNa_M2']
        
        site_occupancy[sample] = min_site

    return site_occupancy       

def initial_parameters(input_csv, ferric_choice):
    '''
    Input your .csv file of xenolith data from the provided template, and this function will output a dictionary containing wfu, so, and bp.

    After invoking function, you should pull out wfu, so, and bp as follows:
    [start code]

    parameters_list = tb_param.initial_parameters('xenolith_small.csv')
    wfu = parameters_list['wfu']
    so = parameters_list['so']
    bp = parameters_list['bp']
    '''
    parameters_list = {}
    
    xeno = pd.read_csv(input_csv)

    working_form_units = formula_units_fast(xeno)

# These are the 3 dictionaries that contain factors for P-T calculations

    wfu = working_form_units
    so = site_occupancy(wfu, ferric_choice)
    bp = beyer_paramaters(wfu, so)


    parameters_list['wfu'] = wfu
    parameters_list['so'] = so
    parameters_list['bp'] = bp

    return parameters_list

def new_template():
    '''
Makes a new xenolith data template .csv file for you to input your xenolith data.
    '''

    if os.path.isfile('xenolith_dataTemplate.csv'):
        print('Data template file already exists in the current directory.')
    else:
        template = [{'Locality': 'Test 1', 'Reference': 'test study', 'SampleID': 'test1', 'cpx_SiO2': '50.58964', 'cpx_TiO2': '0.588516', 'cpx_Al2O3': '8.649199', 'cpx_Cr2O3': '0.07089805', 'cpx_FeO': '3.376068', 'cpx_MnO': '0.03865455', 'cpx_NiO': '0', 'cpx_MgO': '13.091855', 'cpx_CaO': '19.602325', 'cpx_Na2O': '2.4164125', 'cpx_K2O': '0.00064325', 'cpx_Total': '98.42', 'opx_SiO2': '53.11461304', 'opx_TiO2': '0.120747652', 'opx_Al2O3': '5.781704348', 'opx_Cr2O3': '0.038188087', 'opx_FeO': '7.601890435', 'opx_MnO': '0.063495348', 'opx_NiO': '0', 'opx_MgO': '31.73361739', 'opx_CaO': '0.252286652', 'opx_Na2O': '0.041200913', 'opx_K2O': '6.15652E-05', 'opx_Total': '98.75', 'gt_SiO2': '41.11361667', 'gt_TiO2': '0.090340729', 'gt_Al2O3': '24.25045208', 'gt_Cr2O3': '0.058934667', 'gt_FeO': '10.53211667', 'gt_MnO': '0.224152396', 'gt_NiO': '0', 'gt_MgO': '19.20767083', 'gt_CaO': '4.561709375', 'gt_Na2O': '0.02429225', 'gt_K2O': '0', 'gt_Total': '100.06', 'olv_SiO2': '0', 'olv_TiO2': '0', 'olv_Al2O3': '0', 'olv_Cr2O3': '0', 'olv_FeO': '0', 'olv_MnO': '0', 'olv_NiO': '0', 'olv_MgO': '0', 'olv_CaO': '0', 'olv_Na2O': '0', 'olv_K2O': '0', 'olv_Total': '0'},
                    {'Locality': 'Test 2', 'Reference': 'test study', 'SampleID': 'test2', 'cpx_SiO2': '50.93', 'cpx_TiO2': '0.68', 'cpx_Al2O3': '8.86', 'cpx_Cr2O3': '0.02', 'cpx_FeO': '5.45', 'cpx_MnO': '0.033', 'cpx_NiO': '0', 'cpx_MgO': '12.58', 'cpx_CaO': '17.64', 'cpx_Na2O': '3.02', 'cpx_K2O': '0', 'cpx_Total': '99.21', 'opx_SiO2': '52.92', 'opx_TiO2': '0.16', 'opx_Al2O3': '5.14', 'opx_Cr2O3': '0.01', 'opx_FeO': '11.45', 'opx_MnO': '0.09', 'opx_NiO': '0', 'opx_MgO': '29.13', 'opx_CaO': '0.55', 'opx_Na2O': '0.11', 'opx_K2O': '0', 'opx_Total': '99.56', 'gt_SiO2': '40.49', 'gt_TiO2': '0.14', 'gt_Al2O3': '23.85', 'gt_Cr2O3': '0.01', 'gt_FeO': '15.25', 'gt_MnO': '0.33', 'gt_NiO': '0', 'gt_MgO': '16.28', 'gt_CaO': '4.37', 'gt_Na2O': '0.03', 'gt_K2O': '0', 'gt_Total': '100.75', 'olv_SiO2': '0', 'olv_TiO2': '0', 'olv_Al2O3': '0', 'olv_Cr2O3': '0', 'olv_FeO': '0', 'olv_MnO': '0', 'olv_NiO': '0', 'olv_MgO': '0', 'olv_CaO': '0', 'olv_Na2O': '0', 'olv_K2O': '0', 'olv_Total': '0'}]
        dff = pd.DataFrame(template)
        dff.to_csv('xenolith_dataTemplate.csv', index = False)
        print('New template title "xenolith_dataTemplate.csv" should be in current directory!')
    return None

#RUNNING CODE 
