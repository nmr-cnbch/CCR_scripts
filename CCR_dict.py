import sys


CCR_dict = {
    "CCR_1":{
        "ratename":"Ca1Ha1_N2Hn2",
        "TC": 0.0286,
        "angle": "psi",
        "angle_pos": -1,
        "angle_num": 1,
    },
    "CCR_2":{
        "ratename":"N1Hn1_Ca1Ha1",
        "TC": 0.0286,
        "angle":"phi",
        "angle_pos":0,
        "angle_num": 1
        },
    "CCR_3":{
        "ratename":"N1Hn1_N2Hn2",
        "angle":"phi, psi",
        "angle_pos":-1,
        "angle_num": 2
        },
    "CCR_4":{
        "ratename":"Ha1Hn2_cCSA1",
        "angle":"phi",
        "angle_pos":-1,
        "angle_num": 1
        },
    "CCR_5":{
        "ratename":"Ca1Ha1_cCSA1",
        "TC": 0.0286,
        "angle":"psi",
        "angle_pos":-1,
        "angle_num": 1
        },
    "CCR_6":{
        "ratename":"cCSA0_Ca1Ha1",
        "TC": 0.0286,
        "angle":"phi",
        "angle_pos":0,
        "angle_num": 1
        },
    "CCR_7":{
        "ratename":"N1Hn1_cCSA1",
        "angle":"phi, psi",
        "angle_pos":0,
        "angle_num": 2
        },
    "CCR_8":{
        "ratename":"cCSA0_cCSA1",
        "angle":"phi, psi",
        "angle_pos":0,
        "angle_num": 1
        },
    "CCR_9":{
        "ratename":"Ca1Ha1_Ca2Ha2",
        "TC": 0.0286,
        "angle":"phi, psi",
        "angle_pos":-1,
        "angle_num": 2
        }
}


def CCRname2Ratename(ccrname):
    if ccrname=="CCR_1": i="Ca1Ha1_N2Hn2"
    elif ccrname=="CCR_2": i="N1Hn1_Ca1Ha1"
    elif ccrname=="CCR_3": i="N1Hn1_N2Hn2"
    elif ccrname=="CCR_4": i="Ha1Hn2_cCSA1"
    elif ccrname=="CCR_5": i="Ca1Ha1_cCSA1"
    elif ccrname=="CCR_6": i="cCSA0_Ca1Ha1"
    elif ccrname=="CCR_7": i="N1Hn1_cCSA1"
    elif ccrname=="CCR_8": i="cCSA0_cCSA1"
    elif ccrname=="CCR_9": i="Ca1Ha1_Ca2Ha2"
    
    else:
        print ("CCRname2Ratename - Wrong name of CCR (",ccrname,"), check the experiment set file")
        sys.exit(1)
    return i

def CCRname2Ratename2(ccrname):
    if ccrname=="CCR_1": i="Ca1Ha1_N2H2"
    elif ccrname=="CCR_2": i="Ca1Ha1_N1H1"
    elif ccrname=="CCR_3": i="N1H1_N2H2"
    elif ccrname=="CCR_4": i="Ha1H2_C1"
    elif ccrname=="CCR_5": i="Ca1Ha1_C1"
    elif ccrname=="CCR_6": i="Ca1Ha1_C0"
    elif ccrname=="CCR_7": i="N1H1_C1"
    elif ccrname=="CCR_8": i="C0_C1"
    elif ccrname=="CCR_9": i="Ca1Ha1_Ca2Ha2"
    else:
        print ("CCRname2Ratename2 - Wrong name of CCR (",ccrname,"), check the experiment set file")
        sys.exit(1)
    return i


def Ratename2CCRname(ratename):
    if ratename=="Ca1Ha1_N2Hn2": i="CCR_1"
    elif ratename=="N1Hn1_Ca1Ha1": i="CCR_2"
    elif ratename=="N1Hn1_N2Hn2": i="CCR_3"
    elif ratename=="Ha1Hn2_cCSA1": i="CCR_4"
    elif ratename=="Ca1Ha1_cCSA1": i="CCR_5"
    elif ratename=="cCSA0_Ca1Ha1": i="CCR_6"
    elif ratename=="N1Hn1_cCSA1": i="CCR_7"
    elif ratename=="cCSA0_cCSA1": i="CCR_8"
    elif ratename=="Ca1Ha1_Ca2Ha2": i="CCR_9"
    else:
        print ("Ratename2CCRname - Wrong name of CCR rate (",ratename,"), check the experiment set file")
        sys.exit(1)
    return i