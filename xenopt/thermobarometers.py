'''
Mathematical functions for calculating pressure and temperature given the mineral chemistries of garnet, clinopyroxene, and orthopyroxene.

Made by Khalil Droubi, January 2021

Needs to be used with the tb_parameter module.
Before using any of these functions, type:

[start of code]
import tb_parameters as tb_param
import thermobarometers as tb

parameters_list = tb_param.initial_parameters({Your csv file here})

wfu = parameters_list['wfu'] // Weighted formula units for minerals
so = parameters_list['so']   // Site occupancies
bp = parameters_list['bp']   // Separately calculated parameters for the Beyer barometer
[end of code]

Ignore the warnings about double scalars for now.
'''

#imports
import math
import pandas as pd
#import os
#import csv
#import copy
#import matplotlib

# Thermometers:
# Wells 1977 (OLD)
# EG 1979 (OLD)
# Krogh 88 (OLD)

# EG79 w/Fe3+ (PASSED) gt-cpx
# K88 w/Fe3+ (PASSED) gt-cpx
# TBKN 1989 (PASSED) cpx-opx
# T Ca in Opx (PASSED) opx
# Nimis & Taylor (May fail on occasion- math domain error related to certain samples) cpx
# Nakamara(2009) (kind of PASSED? Need to test further) gt-cpx

# Krogh-Ravna(2000) (maybe in future?)
# Putirka (2008) (maybe in future?)

# Barometers
# NG85 (PASSED) gt-opx?
# PBKN (PASSED) gt-opx
# Beyer(2016) (kind of PASSED? Need to test further) gt-cpx

#thermometers

def tbkn(sample_wfu, sample_so, pressure):
    '''
    Clinopyroxene-orthopyroxene Mg-Fe exchange thermometer from Brey and Kohler, 1990

    Brey, G. P., & Köhler, T. (1990). Geothermobarometry in four-phase lherzolites II. New
thermobarometers, and practical assessment of existing thermobarometers. Journal of Petrology,
31(6), 1353-1378.
    '''
    
    # Input T = degrees celsius and P = kilobars
    p = pressure
    
    t25 = (1 - (sample_wfu['cpx_Ca'] / (1 - sample_wfu['cpx_Na']))) / (1 - (sample_wfu['opx_Ca'] / (1 - sample_wfu['opx_Na'])))
    
    numerator = 23664 + (p * (24.9 + (126.3 * sample_so['cpx_Fe#'])))
    denominator = 13.38 + (math.log(t25)**2) + (11.59 * sample_so['opx_Fe#'])
    
    temp_k = (numerator / denominator) #in kelvin
    
    temp_c = temp_k - 273.15 # in celsius
    
    #print ("t25 = " + str(t25))
    
    return temp_c


def T_Ca_in_Opx(sample_wfu, sample_so, pressure):
    '''
    Ca-in-orthopyroxene thermometer from Brey and Kohler, 1990

    Brey, G. P., & Köhler, T. (1990). Geothermobarometry in four-phase lherzolites II. New
thermobarometers, and practical assessment of existing thermobarometers. Journal of Petrology,
31(6), 1353-1378.
    '''
     
    # Input T = degrees celsius and P = kilobars
    p = pressure
    
    temp_k = (6425 + 26.4 * p) / (-math.log(sample_wfu['opx_Ca']) + 1.843)   #in kelvin

    temp_c = temp_k - 273.15    #in celsius
    return temp_c    

def T_K88_wFe3(sample_wfu, sample_so, pressure):
    '''
    Garnet-clinopyroxene Mg-Fe exchange thermometer modified from Krogh, 1988

    Modifications made to correct for Fe3+ in clinopyroxene

    Krogh, E. J. (1988). The garnet-clinopyroxene Fe-Mg geothermometer—a reinterpretation of existing
experimental data. Contributions to Mineralogy and Petrology, 99(1), 44-48.
    '''
    # Input T = degrees celsius and P = kilobars
    p = pressure
    
    temp_k = ( (-6173 * (sample_so['gt_XCa']**2)) + (6731 * sample_so['gt_XCa']) + 1879 + (10 * p) ) / (math.log(sample_so['gt_Fe/Mg'] / sample_so["cpx_Fe2+/Mg"]) + 1.393) #in kelvin
    temp_c = temp_k - 273.15    #in celsius
    
    return temp_c  

def T_EG79_wFe3(sample_wfu, sample_so, pressure):
    '''
    Garnet-clinopyroxene Mg-Fe exchange thermometer modified from Ellis and Green, 1979

    Modifications made to correct for Fe3+ in clinopyroxene

    Ellis, D. J., & Green, D. H. (1979). An experimental study of the effect of Ca upon garnet-clinopyroxene
Fe-Mg exchange equilibria. Contributions to Mineralogy and Petrology, 71(1), 13-22.
    '''
    # Input T = degrees celsius and P = kilobars
    p = pressure
    
    temp_k = (3104 * (sample_wfu['gt_Ca'] / (sample_wfu['gt_Ca'] + sample_wfu['gt_Fe'] + sample_wfu['gt_Mg'])) + 3030 + (10.86 * p) ) / (math.log( (sample_wfu['gt_Fe'] / sample_wfu['gt_Mg']) / ( (sample_wfu['cpx_Fe'] - sample_so['cpx_Fe3+']) / sample_wfu['cpx_Mg'])) + 1.9034) #in kelvin
    
    temp_c = temp_k - 273.15    #in celsius
    
    return temp_c

def T_Nimis_and_Taylor(sample_wfu, sample_so, pressure):
    '''
    Clinopyroxene thermometer from Nakamura, 2009

    Taylor, W. R. (1998). An experimental test of some geothermometer and geobaro-meter formulations for
upper mantle peridotites with application to the ther-mobarometry of fertile lherzolite and garnet
websterite. Neues Jahrbuch für Mineralogie-Abhandlungen, 381-408.

    Nimis, P., & Grütter, H. (2010). Internally consistent geothermometers for garnet peridotites and
pyroxenites. Contributions to Mineralogy and Petrology, 159(3), 411-427.
    '''
    # Input T = degrees celsius and P = kilobars
    p = pressure
    
    temp_k = (23116 + 39.28 * p) / (13.25 + 15.35 * sample_wfu['cpx_Ti'] + 4.5 * sample_wfu['cpx_Fe'] - (1.55 * (sample_wfu['cpx_Al'] + sample_wfu['cpx_Cr'] - sample_wfu['cpx_Na'] - sample_wfu['cpx_K'])) + (math.log(sample_so["cpx_a(en)"])**2) ) 
    
    temp_c = temp_k - 273.15    #in celsius
    
    return temp_c

def T_N09(sample_wfu, sample_so, pressure):
    '''
    Garnet-clinopyroxene thermometer from Nakamura, 2009


    Nakamura, D. (2009). A new formulation of garnet–clinopyroxene geothermometer based on
accumulation and statistical analysis of a large experimental data set. Journal of Metamorphic
Geology, 27(7), 495-508.
    '''
    # Input T = degrees celsius and P = kilobars
    p = pressure
    
    Xprp = sample_so["gt_XMg"]
    Xalm = sample_so["gt_XFe"]
    Xsps = sample_so["gt_XMn"]
    Xgrs = sample_so["gt_XCa"]
    cpx_Xmg = sample_wfu['cpx_Mg'] / (sample_wfu['cpx_Al'] + sample_wfu['cpx_Fe'] + sample_wfu['cpx_Mg'])
    cpx_Xfe = sample_wfu['cpx_Fe'] / (sample_wfu['cpx_Al'] + sample_wfu['cpx_Fe'] + sample_wfu['cpx_Mg'])
    A = 0.5 * Xgrs * (Xprp - Xalm - Xsps)
    B = 0.5 * Xgrs * (Xprp - Xalm + Xsps)
    C = 0.5 * (Xgrs + Xsps) * (Xprp - Xalm)
    Kd = sample_so["gt_Fe/Mg"] / sample_so["cpx_Fe/Mg"]
    
    
    v1 = 2784 + 14.52 * p + (2601 + 1.44 * p) * (2 * Xgrs * Xprp - A)
    v2 = (1183 + 6.98 * p) * (Xgrs ** 2 - A)
    v3 = 105 * (2 * Xgrs * Xalm + B)
    v4 = (814.6 + 3.61 * p) * (Xgrs**2 + B)
    v5 = (254.6 + 8.42 * p) * (2 * Xprp * Xalm - Xalm**2 + C)
    v6 = 83.6 * (Xprp**2 - 2 * Xprp * Xalm + C)
    v7 =  1388 * Xsps - 462 * (cpx_Xmg - cpx_Xfe)
    v_numerator = v1 + v2 - v3 + v4 - v5 - v6 + v7
    u1 = math.log(Kd) + 1.431 + 0.695 * (2 * Xgrs * Xprp + Xgrs**2 - 2 * A)
    u2 = 0.203 * (Xgrs**2 - 2 * Xgrs * Xalm) + 0.922 * Xsps
    
    temp_k = v_numerator / (u1 + u2) # in kelvin
    temp_c = temp_k - 273.15    #in celsius
    
    return temp_c
                                    
#barometers                                    
                                    
def pbkn(sample_wfu, sample_so, sample_bp, temperature):
    '''
    Garnet-orthopyroxene baromometer from Brey and Kohler, 1990

    Brey, G. P., & Köhler, T. (1990). Geothermobarometry in four-phase lherzolites II. New
thermobarometers, and practical assessment of existing thermobarometers. Journal of Petrology,
31(6), 1353-1378.
    '''
    # Input T = degrees celsius and P = kilobars
    
    t = temperature + 273 # in kelvin
    
    Kd = ( ((1 - sample_so['gt_XCa']) ** 3) * (sample_so['gt_XAl'] ** 2) ) / (sample_so['opx_Xmf_M1'] * (sample_so['opx_Xmf_M2'] ** 2) * sample_so['opx_XAl_M1'])
    #Kd PASSED
    
    c1 = ((-8.3143) * t * math.log(Kd)) - 5510 + (88.91 * t) - (19 * (t**1.2)) + (3 * (sample_so['gt_XCa']**2) * 82458) + (sample_so['opx_Xmg_M1'] * sample_so['opx_Xfe_M1'] * (80942 - (46.7 * t))) - (3 * sample_so['gt_XFe'] * sample_so['gt_XCa'] * 17793) - (sample_so['gt_XCa'] * sample_so['gt_XCr'] * ((1.164 * 10**6) - (420.4 * t))) - (sample_so['gt_XFe'] * sample_so['gt_XCr'] * (((-1.25) * 10**6) + (565 * t)))      
    #c1 PASSED
    
    c2 = -0.832 - (8.78 * 10**(-5) * (t - 298)) + (3 * sample_so['gt_XCa']**2 * 3.305) - (sample_so['gt_XCa'] * sample_so['gt_XCr'] * 13.45) + (sample_so['gt_XFe'] * sample_so['gt_XCr'] * 10.5)
    #c2 PASSED
    
    pressure = (-1*c2 - math.sqrt(c2**2 + (4 * 16.6 * 10**(-4) * c1 / 1000))) / (2 * 16.6 * 10**(-4))
    #slightly off by 0.005 Kbar
    
    
    #print("Kd = " + str(Kd))
    #print('c1 = ' + str(c1))
    #print('c2 = ' + str(c2))
    #print(pressure)
    
    return pressure

def NG85(sample_wfu, sample_so, sample_bp, temperature):
    '''
    Garnet-orthopyroxene baromometer from Nickel and Green, 1985


    Nickel, K. G., & Green, D. H. (1985). Empirical geothermobarometry for garnet peridotites and
implications for the nature of the lithosphere, kimberlites and diamonds. Earth and Planetary
Science Letters, 73(1), 158-170.
    '''
    # Input T = degrees celsius and P = kilobars
    
    t = temperature + 273 # in kelvin
    
    NG1 = (-1) / (183.3 + 178.98 * sample_so['opx_XM1_Al_NG'] * (1 - sample_so['opx_XM1_Al_NG']))
    #NG1 PASSED
    
    NG2 = (-1.98717 * t) * math.log( ( ((1 - sample_so['gt_XCa']) ** 3) * (sample_so['gt_XAl'] ** 2) ) / (sample_so['opx_XM1_Fe,Mg'] * (sample_so['opx_Mg+Fe_M2']**2) * sample_so['opx_XM1_Al_NG']) ) - (9000 * (sample_so['gt_XCa']**2)) - 3400 * (2 * (sample_so['gt_XCr']**2) - (sample_so['opx_XM1_Mg'] * sample_so['opx_XM1_Cr'])) - (sample_so['gt_XCa'] * sample_so['gt_XCr'] * (90853 - (52.1 * t))) - (7590 * sample_so['gt_XFe'] * sample_so['gt_XCa']) + (5157 * sample_so['opx_XM1_Mg'] * sample_so['opx_XM1_Fe']) + 6047 - (3.23 * t)      
    #NG2 PASSED
    
    pressure = NG1 * NG2
    
    #print(NG1)
    #print(NG2)
    
    return pressure


def P_Beyer(sample_wfu, sample_so, sample_bp, temperature):
    '''
    Garnet-clinopyroxene "eclogite" barometer from Beyer et al., 2015

    Beyer, C., Frost, D. J., & Miyajima, N. (2015). Experimental calibration of a garnet–clinopyroxene
geobarometer for mantle eclogites. Contributions to Mineralogy and Petrology, 169(2), 1-21.
    '''


    # Input T = degrees celsius and P = kilobars
    
    t = temperature + 273 # in kelvin
    
 # Interaction Parameters
    #garnet Margules parameters, values from Ganguly et al., 1996, with pressure dependent terms equal to zero
    gt_W_CaMg = 21627 - (t * 5.78)
    gt_W_MgCa = 9834 - (t * 5.78)
    gt_W_CaFe = 873 - (1.69 * t)
    gt_W_FeCa = 6773 - (1.69 * t)
    gt_W_MgFe = 2117
    gt_W_FeMg = 695

    #gt_grs_dodec_RTg = ((0.5 * gt_W_CaMg * sample_bp['gt_XMg'] * (1 - sample_bp['gt_XCa'] + sample_bp['gt_XMg'] + 2 * sample_bp['gt_XCa'] * (sample_bp['gt_XCa'] - sample_bp['gt_XMg'] - 1)) )  + (0.5 * gt_W_MgCa * sample_bp['gt_XMg'] * ( 1 - sample_bp['gt_XCa'] - sample_bp['gt_XMg'] - (2 * sample_bp['gt_XCa'] * (sample_bp['gt_XCa'] - sample_bp['gt_XMg'] - 1)) + 0.5 * gt_W_CaFe * sample_bp['gt_XFe'] * (1 - sample_bp['gt_XCa'] + sample_bp['gt_XFe'] + (2 * sample_bp['gt_XCa'] * (sample_bp['gt_XCa'] - sample_bp['gt_XFe'] -1))) + (0.5 * gt_W_FeCa * sample_bp['gt_XFe'] * (1 - sample_bp['gt_XCa'] - sample_bp['gt_XFe'] - (2 * sample_bp['gt_XCa'] * (sample_bp['gt_XCa'] - sample_bp['gt_XFe'] - 1)))) + (gt_W_MgFe * sample_bp['gt_XMg'] * sample_bp['gt_XFe'] * (sample_bp['gt_XFe'] - sample_bp['gt_XMg'] - 0.5)) + (gt_W_FeMg * sample_bp['gt_XFe'] * sample_bp['gt_XMg'] * (sample_bp['gt_XFe'] - sample_bp['gt_XMg'] - 0.5)) + (  (((gt_W_MgCa - gt_W_CaMg) + (gt_W_CaFe - gt_W_FeCa) + (gt_W_FeMg - gt_W_MgFe)) / 2) * sample_bp['gt_XMg'] * sample_bp['gt_XFe'] * (1 - 2 * sample_bp['gt_XCa']) ) )) ) 
    
    p1 = 0.5 * gt_W_CaMg * sample_bp['gt_XMg'] * (1 - sample_bp['gt_XCa'] + sample_bp['gt_XMg'] + 2 * sample_bp['gt_XCa'] * (sample_bp['gt_XCa'] - sample_bp['gt_XMg'] - 1))             
    p2 = 0.5 * gt_W_MgCa * sample_bp['gt_XMg'] * (1 - sample_bp['gt_XCa'] - sample_bp['gt_XMg'] - 2 * sample_bp['gt_XCa'] * (sample_bp['gt_XCa'] - sample_bp['gt_XMg'] - 1))                                                                                                                               
    p3 = 0.5 * gt_W_CaFe * sample_bp['gt_XFe'] * (1 - sample_bp['gt_XCa'] + sample_bp['gt_XFe'] + 2 * sample_bp['gt_XCa'] * (sample_bp['gt_XCa'] - sample_bp['gt_XFe'] -1))                                        
    p4 = 0.5 * gt_W_FeCa * sample_bp['gt_XFe'] * (1 - sample_bp['gt_XCa'] - sample_bp['gt_XFe'] - 2 * sample_bp['gt_XCa'] * (sample_bp['gt_XCa'] - sample_bp['gt_XFe'] - 1))                                        
    p5 =       gt_W_MgFe * sample_bp['gt_XMg'] * sample_bp['gt_XFe'] * (sample_bp['gt_XFe'] - sample_bp['gt_XMg'] - 0.5)                                        
    p6 =       gt_W_FeMg * sample_bp['gt_XFe'] * sample_bp['gt_XMg'] * (sample_bp['gt_XFe'] - sample_bp['gt_XMg'] - 0.5)
    p7 =  (((gt_W_MgCa - gt_W_CaMg) + (gt_W_CaFe - gt_W_FeCa) + (gt_W_FeMg - gt_W_MgFe)) / 2) * sample_bp['gt_XMg'] * sample_bp['gt_XFe'] * (1 - 2 * sample_bp['gt_XCa'])   
    
    gt_grs_dodec_RTg = p1 + p2 + p3 + p4 + p5 + p6 + p7
    gt_grs_RTLN_X3 = 8.314 * t * math.log(sample_bp['gt_XCa']**3)
    gt_grs_Rtlna_CaAl = 3 * gt_grs_dodec_RTg + gt_grs_RTLN_X3
    
    u1 = 0.5 * gt_W_MgCa * sample_bp['gt_XCa'] * (1 - sample_bp['gt_XMg'] + sample_bp['gt_XCa'] + 2 * sample_bp['gt_XMg'] * (sample_bp['gt_XMg'] - sample_bp['gt_XCa'] - 1))             
    u2 = 0.5 * gt_W_CaMg * sample_bp['gt_XCa'] * (1 - sample_bp['gt_XMg'] - sample_bp['gt_XCa'] - 2 * sample_bp['gt_XMg'] * (sample_bp['gt_XMg'] - sample_bp['gt_XCa'] - 1))                                                                                                                               
    u3 = 0.5 * gt_W_MgFe * sample_bp['gt_XFe'] * (1 - sample_bp['gt_XMg'] + sample_bp['gt_XFe'] + 2 * sample_bp['gt_XMg'] * (sample_bp['gt_XMg'] - sample_bp['gt_XFe'] -1))                                        
    u4 = 0.5 * gt_W_FeMg * sample_bp['gt_XFe'] * (1 - sample_bp['gt_XMg'] - sample_bp['gt_XFe'] - 2 * sample_bp['gt_XMg'] * (sample_bp['gt_XMg'] - sample_bp['gt_XFe'] - 1))                                        
    u5 =       gt_W_CaFe * sample_bp['gt_XCa'] * sample_bp['gt_XFe'] * (sample_bp['gt_XCa'] - sample_bp['gt_XFe'] - 0.5)                                        
    u6 =       gt_W_FeCa * sample_bp['gt_XFe'] * sample_bp['gt_XCa'] * (sample_bp['gt_XFe'] - sample_bp['gt_XCa'] - 0.5)
    u7 =  (((gt_W_MgCa - gt_W_CaMg) + (gt_W_CaFe - gt_W_FeCa) + (gt_W_FeMg - gt_W_MgFe)) / 2) * sample_bp['gt_XCa'] * sample_bp['gt_XFe'] * (1 - 2 * sample_bp['gt_XMg'])   
    
    gt_pyp_dodec_RTg = u1 + u2 + u3 + u4 + u5 + u6 + u7
    gt_pyp_RTLN_X3Y2 = 8.314 * t * math.log(sample_bp['gt_XMg']**3)
    gt_pyp_Rtlna_MgAl = 3 * gt_pyp_dodec_RTg + gt_pyp_RTLN_X3Y2
    
        #clinopyroxene Margules parameters, values from Ganguly et al., 1996
    
    cpx_W_MgFe_M1 = 0
    cpx_W_MgAl_M1 = 4849
    cpx_W_FeAl_M1 = 3573
    cpx_W_CaNa_M2 = 11567
    cpx_W_MgNa_M2 = -35873
    cpx_W_FeNa_M2 = 59319
    cpx_W_CaMg_M2 = 24920
    cpx_W_CaFe_M2 = -49275
    cpx_W_FeMg_M2 = 0
    cpx_W_AlSi_M2 = 7038
    
    cpx_diop_RTln = 8.314 * t * math.log(sample_bp['cpx_XMg_M1'] * sample_bp['cpx_XCa_M2'] * (sample_wfu['cpx_Si'] / 2)**2 )
    
    v1 = (sample_bp['cpx_XMg_M2'] ** 2) * cpx_W_CaMg_M2 + (sample_bp['cpx_XNa_M2'] ** 2) * cpx_W_CaNa_M2 + (sample_bp['cpx_XFe_M2'] ** 2) * cpx_W_CaFe_M2
    v2 = sample_bp['cpx_XMg_M2'] * sample_bp['cpx_XFe_M2'] * (cpx_W_CaMg_M2 + cpx_W_CaFe_M2 -  cpx_W_FeMg_M2) + sample_bp['cpx_XMg_M2'] * sample_bp['cpx_XNa_M2'] * (cpx_W_CaMg_M2 + cpx_W_CaNa_M2 - cpx_W_MgNa_M2) + sample_bp['cpx_XNa_M2'] * sample_bp['cpx_XFe_M2'] * (cpx_W_CaNa_M2 + cpx_W_CaFe_M2 -  cpx_W_FeNa_M2)
    cpx_diop_RTlna_M2 = v1 + v2
    
    cpx_diop_RTlna_M1 = sample_bp['cpx_XAl_M1']**2 * cpx_W_MgAl_M1 + sample_bp['cpx_XFe3+_M1']**2 * cpx_W_MgFe_M1 + sample_bp['cpx_XAl_M1'] * sample_bp['cpx_XFe3+_M1'] * (cpx_W_MgAl_M1 + cpx_W_MgFe_M1 - cpx_W_FeAl_M1)
    
    cpx_diop_RTlna_TET = cpx_W_AlSi_M2 * ((sample_bp['cpx_Al_IV'] / 2)**2)
    
    cpx_diop_reciprocal = (-sample_bp['cpx_XAl_M1']) * sample_bp['cpx_XFe_M2'] * (-21000) - (sample_bp['cpx_XFe3+_M1'] * sample_bp['cpx_XFe_M2'] * (-4680)) + (sample_bp['cpx_XMg_M2'] * (-5010))
    # FAILED?
    cpx_diop_RTlna_CaMgSi = cpx_diop_RTln + cpx_diop_RTlna_M2 + cpx_diop_RTlna_M1 + cpx_diop_RTlna_TET + cpx_diop_reciprocal
                                                                                    
    cpx_CaTs_RTln = 8.314 * t * math.log(sample_bp['cpx_XAl_M1'] * (4 * sample_bp['cpx_XCa_M2']) * (sample_bp['cpx_Al_IV'] / 2) * (sample_wfu['cpx_Si'] / 2) )
    
    cpx_CaTs_RTlna_M2 = cpx_diop_RTlna_M2
    cpx_CaTs_RTlna_M1 = sample_bp['cpx_XMg_M1']**2 * cpx_W_MgAl_M1 + sample_bp['cpx_XFe3+_M1']**2 * cpx_W_FeAl_M1 + sample_bp['cpx_XMg_M1'] * sample_bp['cpx_XFe3+_M1'] * (cpx_W_MgAl_M1 +  cpx_W_FeAl_M1 - cpx_W_MgFe_M1)
    cpx_CaTs_RTlna_TET = cpx_W_AlSi_M2 * (sample_wfu['cpx_Si'] / 2)**2
    
    cpx_CaTs_reciprocal =  (1 - sample_bp['cpx_XAl_M1']) * sample_bp['cpx_XFe_M2'] * (-21000) - (sample_bp['cpx_XFe3+_M1'] * sample_bp['cpx_XFe_M2'] * (-4680)) + (sample_bp['cpx_XMg_M2'] * (-5010))
    # FAILED?
    
    cpx_CaTs_RTlna_CaAlAl = cpx_CaTs_RTln + cpx_CaTs_RTlna_M2 + cpx_CaTs_RTlna_M1 + cpx_CaTs_RTlna_TET + cpx_CaTs_reciprocal
    
#     print(cpx_diop_reciprocal)
#     print(cpx_CaTs_reciprocal)
#     print(cpx_CaTs_RTlna_CaAlAl)
    
    Eq_constant_RTlnK = cpx_CaTs_RTlna_CaAlAl + cpx_diop_RTlna_CaMgSi - ((1/3) * gt_pyp_Rtlna_MgAl) - ((2/3) *  gt_grs_Rtlna_CaAl)
    #About 800 off, but is this due to calculation error or differences in Fe3+?
    
    Eq_constant_K = math.exp(Eq_constant_RTlnK / (8.314 * t))
    # off by 0.16
   
    #Fitting parameters from Beyer (2015)
    delta_S = 23.5
    delta_H = 16500
    delta_Vr = 0.719
    
    delta_Gt = delta_H - delta_S * t
    
    pressure = (-(delta_Gt + 8.314 * t * math.log(Eq_constant_K))) / (delta_Vr * 1000)
    
    
    return pressure
    
