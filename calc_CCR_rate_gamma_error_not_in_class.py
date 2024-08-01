#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 27 11:18 2022

@author: Paulina Bartosinska-Marzec
"""
""" The script to reading Sparky peak list and calculating CCR rates"""
    
    
import array
from genericpath import exists
from math import atanh
from copy import deepcopy
import math
import sys
import csv
from typing import Union
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy import stats
from sklearn.linear_model import LinearRegression
import os
import logging

if len(sys.argv)==1:
    sys.exit("""
    The script to calculate CCR rates. In director has to contains:
    - peak lists (with peak position in points and peak height)
    - experiment_set
    - sequence in FASTA format

    For run script type in command line:
        python3 calc_CCR_rate.py [file director]

    additionally you can add:
        --seq  - if name of file with amino acid sequence is not 'seq', add this with name of file
        --refgamma  - if you have file with reference values of CCR rates, add this with name of file
                      (file must be .csv, columns name should be: AA, psi_angle (or/and phi_angle), CCR name
        --expset    - if you want use experiments setup file with different filename than "experiments_set.txt" (structure of file must be the same as orginal file) 
        --peak_uncert    - if you want to use peaks uncertainty to calculate CCR rates uncertainity insted of global averenge noise level seted up in experiments_set.txt  """)


"""Command line reading"""

file_director = sys.argv[1]
if file_director[-1] != "/":
    file_director += "/"
    print ("\nFile director",file_director)

if "--seq" in sys.argv:    
    i = sys.argv.index("--seq")
    seq_file_name = sys.argv[i+1] 
else: seq_file_name = "{}seq".format(file_director)
print("File with amino acid sequence: {}".format(seq_file_name)) 

if "--refgamma" in sys.argv: 
    i = sys.argv.index("--refgamma")
    refgammaFlag = True
    gamma_cal_file_name = sys.argv[i+1]
    print ("Reference CCR rates are inclued from: {}".format(gamma_cal_file_name)) 
else: 
    refgammaFlag = False
    gamma_cal_file_name = None

RaportDir = file_director+"all_outputs/"
RaportBoxFile = RaportDir+"RaportBox.txt"
if not os.path.exists(RaportDir):
        os.mkdir(RaportDir)


if "--expset" in sys.argv:    
    i = sys.argv.index("--expset")
    exp_file_name = "{}{}".format(file_director,sys.argv[i+1])
    print("Diffrent experiment set file: {}".format(exp_file_name)) 
else: exp_file_name = "{}experiments_set.txt".format(file_director)

if "--peak_uncert" in sys.argv: 
    PeakUncertiFlag = True
else:
    PeakUncertiFlag = False
# 
if "--publication" in sys.argv: 
    PublicationFlag = True
else:
    PublicationFlag = False

aminoacids=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    
CCR_dict = {"CCR_1":{"angle":"psi","angle_pos":-1,},
            "CCR_2":{"angle":"phi","angle_pos":0,},
            "CCR_3":{"angle":"phi, psi","angle_pos":-1,},
            "CCR_4":{"angle":"phi","angle_pos":-1,},
            "CCR_5":{"angle":"psi","angle_pos":-1,},
            "CCR_6":{"angle":"phi","angle_pos":0,},
            "CCR_7":{"angle":"phi, psi","angle_pos":0,},
            "CCR_8":{"angle":"phi, psi","angle_pos":0,},
            "CCR_9":{"angle":"phi, psi","angle_pos":-1,},
            }

"""Classes"""

class CSpectrum:
    def __init__(self):
        self.CCR_name = ""       # name of CCR: CCR_1, CCR_2,                           -> required
        self.auto_name = ""       # name of file with auto version                      -> required
        self.cross_name = ""       # name of file with cross version                    -> required
        self.n_dim = 0            # number of dimensions                                -> required
        self.nucl_name = []            # nuclei of all peaks: H, N, C, CA, CB, HA, HB
        self.nucl_pos = []             # nuclei position of all peaks: -2, -1, 0, 1, 2 
        
        self.n_angle = 0          # number of measured angle  
        self.angle = []           # angle: phi, psi                                     
        self.angle_pos = []       # angle position of all peaks: -2, -1, 0, 1, 2        
        self.CCR_pos = 100
        self.ns = [0,0]          # number of scans in auto [0] and cross [1] version    -> required
        self.tc_vol = -1.0        # time of ccr evolution                               -> required
        self.Hroi = [0.0,0.0]       # downfield and upfield of direct dimension
        self.noise = [0.0,0.0]       # noise average level auto [0] and cross [1] version    -> required

        self.peak = [0.0,0.0]            #information about peaks from CResidue class
        self.is_peaklist = False
        self.peak_uncertainty = False
        self.other = ""
    
                    
class CResidue:
    def __init__(self):
        self.aa_name = "None"
        self.aa_number = 0              # the number of the amino acid to which this peak corresponds
        self.is_peak = False       # information about the presence of a peak


        self.peak_pos = []            # chemical shifts for all nuclei of peak, length depends on dimentionality
        self.peak_pos_points = []
        self.peak_intens = [0.0,0.0]      # peak hight in auto [0] and cross [1] version
        self.descript = ""              # description of peak which is in first column in Sparky-like peak list
        # self.peak_seq_pos_a = []
        self.is_overlap = False       # information about the presents of overlapping of a peak: maybe, yes, no   
        self.overlap_peaks = {}       #  
        self.peak_uncertainty = [1.0,1.0] 

        self.is_ccr_rate = False
        self.ccr_rate = -99999.0
        self.ccrrate_error_value = -99999.0
        self.ccrrate_calculation_error = False

        self.is_gamma_calc = False      # information about the presence of a reference gamma
        self.gamma_ref = -99999.0       # value of gamma (CCR rate) which will use as reference data for experimental data, value additional file with calculated previosly gammas
        self.Ix_theor = 0.0
        self.Ix_theor_Ia_ratio = 0.0
        self.Ix_theor_Ia_ratio_without_NS = 0.0

        
        

"""Small functions"""


def Res3to1(res):
    if res=="PRO": i="P"
    elif res=="ALA": i="A"
    elif res=="SER": i="S"
    elif res=="THR": i="T"
    elif res=="CYS": i="C"
    elif res=="CYSox": i="C"
    elif res=="CYSred": i="C"
    elif res=="GLY": i="G"
    
    elif res=="ASN": i="N"
    elif res=="ASP": i="D"
    elif res=="PHE": i="F"
    elif res=="TYR": i="Y"
    
    elif res=="ILE": i="I"
    elif res=="LEU": i="L"
    
    elif res=="ARG": i="R"
    elif res=="GLN": i="Q"
    elif res=="GLU": i="E"
    elif res=="VAL": i="V"
    elif res=="LYS": i="K"
    elif res=="MET": i="M"
    
    elif res=="TRP": i="W"
    elif res=="HIS": i="H"
    
    else:
        print ("Wrong name of residue (",res,"), check the seq file - FASTA format")
        sys.exit(1)
    return i

def Res1to3(res):
    if res=="P": i="PRO"
    elif res=="A": i="ALA"
    elif res=="S": i="SER"
    elif res=="T": i="THR"
    elif res=="C": i="CYS"
    elif res=="C_ox": i="CYSox"
    elif res=="C_red": i="CYSred"
    elif res=="G": i="GLY"
    
    elif res=="N": i="ASN"
    elif res=="D": i="ASP"
    elif res=="F": i="PHE"
    elif res=="Y": i="TYR"
    
    elif res=="I": i="ILE"
    elif res=="L": i="LEU"
    
    elif res=="R": i="ARG"
    elif res=="Q": i="GLN"
    elif res=="E": i="GLU"
    elif res=="V": i="VAL"
    elif res=="K": i="LYS"
    elif res=="M": i="MET"
    
    elif res=="W": i="TRP"
    elif res=="H": i="HIS"
    
    else:
        print ("Wrong name of residue (",res,"), check the seq file")
        sys.exit(1)
    return i


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
        print ("Wrong name of CCR (",ccrname,"), check the experiment set file")
        sys.exit(1)
    return i



def extract_number(given_string):
    only_number = ""
    for indexi, i in enumerate(given_string):
        if i.isdigit():
            only_number += i
    return only_number


def changeFlag(flaga):
    if flaga == True:
        flaga = False
    elif flaga == False:
        flaga = True
    return flaga

def calc_distance(peak1,peak2):
    sum_of_squares = 0
    for i in range(len(peak1)):
        # print (peak1, "and", peak2)
        sum_of_squares += (peak1[i]-peak2[i])**2
    dis =  sum_of_squares**(1/len(peak1))
    return dis

def print_raport(text_to_write):
    print (text_to_write)
    RaportBox.write(text_to_write)
    RaportBox.write("\n")



# def OnePeakInfoFromAllExp(peak_num,item_name):
#     list_of_elements = []

#     for experiment in Experiments:
#         list_of_elements.append(deepcopy(experiment.peak[peak_num]."{}".format(item_name)))
#         print ("list_of_elements: peak number = {}, item name = {} ---> {}".format(peak_num, item_name,experiment.peak[peak_num]."{}".format(item_name)))
#     return list_of_elements



"""Reading functions"""


def ReadExpSet(exp_file):               # wczytywanie danych z pliku experiments_set.txt do klasy CSpectrum
    experiments = []
    to_compere_dict = {}
    with open(exp_file, "r") as exp_set:
        # print ("otwarte")
        lines = exp_set.readlines()
        expset_lines=[]
        commentFlag = False
        for indexl, line in enumerate(lines):
            if commentFlag == False:
                if "type_of_CCR" in line:
                    expset_lines.append(deepcopy(indexl))
            if "====" in line:
                expset_lines.append(deepcopy(indexl))
                commentFlag = changeFlag(commentFlag)
        # print ("expset_lines", expset_lines)
        if commentFlag == False:
            expset_lines.append(deepcopy(len(lines)))

        for i in range(0,len(expset_lines)-1):
            # print ("Lines between ", expset_lines[i], expset_lines[i+1])
            one_experinet = CSpectrum()
            # one_experinet.name=lines[expset_lines[i]][8:-1]     # ???????????
            for line in range(expset_lines[i],expset_lines[i+1]):
                if lines[line] != "\n":
                    items=lines[line].split()
                    if "auto_name"in lines[line]:
                        # print ("auto_name_file", items[1])
                        one_experinet.auto_name =items[1]
                    if "cross_name"in lines[line]:
                        # print ("cross_name_file", items[1])
                        one_experinet.cross_name =items[1]
                    if "type_of_CCR" in lines[line]:
                        one_experinet.CCR_name=items[1]
                        if items[1] not in to_compere_dict:
                            to_compere_dict[items[1]]=[i]
                        else:
                            to_compere_dict[items[1]].append(deepcopy(i))
                        # print ("CCR_NAME", items[1])
                    if "dir_auto" in lines[line]:
                        auto_name_dir = items[1]+"_"+one_experinet.CCR_name+"_a"
                        # print ("auto_name_dir", auto_name_dir)
                        one_experinet.auto_name = auto_name_dir
                    if "dir_cross" in lines[line]:
                        cross_name_dir = items[1]+"_"+one_experinet.CCR_name+"_x"
                        # print ("cross_name_dir", cross_name_dir)
                        one_experinet.cross_name=cross_name_dir
                    if "dimension" in lines[line]:
                         one_experinet.n_dim=int(items[1])
                    if "nucl" in lines[line]:
                        for a in range(1, len(items)):
                            if items[a] == "#" : break
                            if items[a] in ["H", "N", "CO", "CA", "CB", "HA", "HB"]:
                                one_experinet.nucl_name.append(deepcopy(items[a]))
                            
                    if "pos_nucl" in lines[line]:
                        for a in range(1, len(items)):
                            if items[a] == "#" : break
                            one_experinet.nucl_pos.append(deepcopy(int(items[a])))
                    if "angle_num" in lines[line]:
                        one_experinet.n_angle=int(items[1])
                    if "angle_pos" in lines[line]:
                        for a in range(1, len(items)):
                            if items[a] == "#" : break
                            one_experinet.angle_pos.append(deepcopy(int(items[a])))
                    if "angle_name" in lines[line]:
                        for a in range(1, len(items)):
                            if items[a] == "#" : break
                            one_experinet.angle.append(deepcopy(items[a]))
                    if "NS_auto" in lines[line]:
                        one_experinet.ns[0]=int(items[1])
                    if "NS_cross" in lines[line]:
                        one_experinet.ns[1]=int(items[1])
                    if "TC" in lines[line]:
                        one_experinet.tc_vol=float(items[1])
                    if "H_roi" in lines[line]:
                        one_experinet.Hroi[0]=float(items[1])
                        one_experinet.Hroi[1]=float(items[2])
                    if "other" in lines[line]:
                        one_experinet.other=items[1]
                        # print ("CCR_NAME", items[1])
                    if "noise" in lines[line]:
                        one_experinet.noise[0]=float(items[1])
                        one_experinet.noise[1]=float(items[2])
            one_experinet.CCR_pos = deepcopy(CCR_dict[one_experinet.CCR_name]["angle_pos"])
            experiments.append(deepcopy(one_experinet))
            print ("{} - {}\n\t reference exp: {}\n\t transfer exp: {}".format(one_experinet.CCR_name,one_experinet.other,one_experinet.auto_name,one_experinet.cross_name))
            # print ("Experiment number:{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}".format(i, one_experinet.CCR_name, one_experinet.n_dim, one_experinet.nucl_name, one_experinet.nucl_pos, one_experinet.n_angle, one_experinet.angle_pos, one_experinet.angle, one_experinet.ns, one_experinet.tc_vol, one_experinet.Hroi), file=RaportBox)
            RaportBox.write(f"""
Experiment number: {i}
CCR name: {one_experinet.CCR_name}
n dim: {one_experinet.n_dim}
angle pos: {one_experinet.angle_pos}
ns: {one_experinet.ns}
Tc: {one_experinet.tc_vol}

""")
    return experiments, to_compere_dict


def ReadSequnce(seq_file_name):
    number_of_aa_per_residue = {}
    sequence = []
    with open(seq_file_name, "r") as input_file:
        linia_seq=input_file.readlines()
        FASTAFlag=False
        FASTAseq = ""
        if ">" in linia_seq[0]:
            FASTAFlag=True
            print ("FASTAFlag")
            # FASTAseq=[]
            aa_no = 0
            for indexl,l_seq in enumerate(linia_seq[1:]):
                FASTAseq += l_seq
                oneletter=list(str(l_seq))
                for aa in oneletter: 
                    if aa.isalpha():
                        seq = CResidue()
                        seq.aa_name = aa         
                        aa_no += 1
                        seq.aa_number = aa_no
                        sequence.append(deepcopy(seq))
                        if aa in number_of_aa_per_residue:
                            number_of_aa_per_residue[aa]+=1
                        else:
                            number_of_aa_per_residue[aa]=1
                        # FASTAseq.append(deepcopy(seq))
            print_raport ("FASTA sequence: {}".format(FASTAseq))

    return sequence,number_of_aa_per_residue

def read_aminoacid_number(peak_description:str):
    aminoacids = peak_description.split("-")
    first_dim_sign_list = list(aminoacids[0])
    last_dim_sign_list = list(aminoacids[-1])
    print_raport(f"{first_dim_sign_list}, {last_dim_sign_list}")
    digit_list_first = []
    digitFlag_first = False
    digit_list_last = []
    digitFlag_last = False

    for indexn, sign in enumerate(first_dim_sign_list):
        if digitFlag_first==False:
            if sign.isdigit():
                digit_list_first.append(deepcopy(indexn))
                digitFlag_first = True
        elif digitFlag_first:
            if not sign.isdigit():
                digit_list_first.append(deepcopy(indexn))
                digitFlag_first = False

    for indexn, sign in enumerate(last_dim_sign_list):
        if digitFlag_last==False:
            if sign.isdigit():
                digit_list_last.append(deepcopy(indexn))
                digitFlag_last = True
        elif digitFlag_last:
            if not sign.isdigit():
                digit_list_last.append(deepcopy(indexn))
                digitFlag_last = False

    if len(digit_list_first)>=2:
        aminoacids_number_first = ''.join([n for n in first_dim_sign_list[digit_list_first[0]:digit_list_first[1]]])
    else: aminoacids_number_first = 'nan'
    
    if len(digit_list_last)>=2:
        aminoacids_number_last = ''.join([n for n in last_dim_sign_list[digit_list_last[0]:digit_list_last[1]]])
    else: aminoacids_number_last = 'nan'
    
    print_raport(f"digit_list_first: {digit_list_first}, digit_list_last: {digit_list_last}")
    print_raport(f"aminoacids_number_first: {aminoacids_number_first}, aminoacids_number_last: {aminoacids_number_last}")

    if aminoacids_number_first == aminoacids_number_last:
        aminoacids_number = aminoacids_number_first
        HaveNumberFlag = True
    elif aminoacids_number_first.isdigit() and aminoacids_number_last == 'nan':
        aminoacids_number = aminoacids_number_first
        HaveNumberFlag = True
    elif aminoacids_number_last.isdigit() and aminoacids_number_first == 'nan':
        aminoacids_number = aminoacids_number_last
        HaveNumberFlag = True
    else:
        if first_dim_sign_list[digit_list_first[1]+1] == 'N' or first_dim_sign_list[digit_list_first[1]+1] == 'H':
            aminoacids_number = aminoacids_number_first
            HaveNumberFlag = True
        elif last_dim_sign_list[digit_list_last[1]+1] == 'N' or last_dim_sign_list[digit_list_last[1]+1] == 'H':
            aminoacids_number = aminoacids_number_last
            HaveNumberFlag = True
        else: HaveNumberFlag = False

    if HaveNumberFlag == False:
        print_raport(f"Probelm with number reading in: {first_dim_sign_list} and {last_dim_sign_list}")
        aminoacids_number = "0"

    print(f"aminoacids_number: {aminoacids_number}")
    return aminoacids_number

def Read_peaklist(file_director, peak_list_name_auto, peak_list_name_cross, s_dim, p_list, is_there_peaklist, points_mode=False):
    
    NameFlag = [False, False]
    auto_noise = 0.0
    cross_noise = 0.0

    # NotFoundName = []
    if points_mode == False:
        listofnames = [".list","_new_ppm.list"]
    else:
        listofnames = ["_points.list","_new_points.list"]
    peak_list_basic_names = [peak_list_name_auto, peak_list_name_cross]
    peak_list_names = ["None","None"]
    for indexlv, list_verson in enumerate(peak_list_basic_names):
        for i, l in enumerate(listofnames):
            peaklistfile = f"{file_director}peak_lists/{list_verson}{l}"
            try:
                f=open(peaklistfile)
                NameFlag[indexlv] = True
                peak_list_names[indexlv] = peaklistfile
                # print ("peaklist --->",peaklistfile, "Done")
                f.close
            except FileNotFoundError:
                # NotFoundName.append(deepcopy(peaklistfile))
                # print ("There is no such file or directory:", peaklistfile, file=RaportBox)
                RaportBox.write("There is no such file or directory: {}\n".format(peaklistfile))
    
    

    if NameFlag[0] == True and NameFlag[1] == True:
        RaportBox.write("\n\nLISTA:{}\n".format(peak_list_names))
        with open(peak_list_names[0], 'r') as pl_a, open(peak_list_names[1], 'r') as pl_x:  
            is_there_peaklist = True
            p_lines_a = pl_a.readlines()
            p_lines_x = pl_x.readlines()
            if points_mode:
            # if points_mode and "Average noise level" in p_lines_a[1] and p_lines_x[1]:
                auto_noise = float(p_lines_a[1].split("=")[1])
                cross_noise = float(p_lines_x[1].split("=")[1])
                print_raport(f"auto_noise = {auto_noise}, cross_noise = {cross_noise}")
                # RaportBox.write(f"auto_noise = {auto_noise}, cross_noise = {cross_noise}")
            for indexl, line in enumerate(p_lines_a):
                if indexl > 1 and len(line) > 1:
                    item_a = line.split()
                    item_x = p_lines_x[indexl].split()
                    # Reading description
                    if item_a[0] != item_x[0]:
                        print ("Error: Auto ({}) and cross ({}) peak list are not compatible:\n{} ? {}".format(peak_list_names[0], peak_list_names[1], item_a[0], item_x[0]))
                        sys.exit()
                    description = item_a[0]
                    aminoacids_number = read_aminoacid_number(description)
                    # aminoacids = description.split("-")
                    # aminoacids_number = ''.join([n for n in str(aminoacids[s_dim-1]) if n.isdigit()])
                    # if aminoacids_number == '':
                    #     aminoacids_number = ''.join([n for n in str(aminoacids[0]) if n.isdigit()])
                    RaportBox.write(f"aminoacids_number: {aminoacids_number,description}\n")
                    aminoacids_number = int(aminoacids_number)
                    for residue in p_list:
                        # residue = CResidue()
                        if aminoacids_number == residue.aa_number:
                            residue.is_peak = True
                            residue.descript = description
                            if points_mode == True:  # for points mode
                                for i in range(1,s_dim+1):
                                    residue.peak_pos_points.append(deepcopy(int(item_a[i])))
                                    # print ("residua append --->", residue.peak_pos_points)
                            else: 
                                for i in range(1,s_dim+1):
                                    residue.peak_pos.append(deepcopy(float(item_a[i])))
                    # Reading peak intensity
                            if len(item_a)>s_dim+1:
                                RaportBox.write("auto = {}, cross = {}\n".format(item_a[s_dim+1],item_x[s_dim+1]))
                                residue.peak_intens[0] = float(item_a[s_dim+1])
                                residue.peak_intens[1] = float(item_x[s_dim+1])
    else:
        print_raport ("Missing file for {} and {}! {}\n".format(peak_list_name_auto, peak_list_name_cross, peak_list_names))
        is_there_peaklist = False
    if points_mode:
        return p_list, is_there_peaklist, auto_noise, cross_noise
    else:
        return p_list, is_there_peaklist

def Read_Reference_Gamma(gamma_cal_file_name):
    with open(gamma_cal_file_name, newline='') as gc_file:  
        reader = csv.DictReader(gc_file, delimiter=",")
        headers = next(reader)
        gamma_file_dict = {}
        for item in headers:
            gamma_file_dict[item]=[headers[item]]
        for col in reader:
            for item in headers:
                # print (item)
                gamma_file_dict[item].append(col[item])
    # print(gamma_file_dict)
    try:
        Res3to1(gamma_file_dict["AA"][0][-3:])
        if gamma_file_dict["AA"][0][:-3].isnumeric():
            gamma_file_dict["seq_num"]=[gamma_file_dict["AA"][0][:-3]]
            for i in range(1,len(gamma_file_dict["AA"])):
                print_raport ("gamma_file_dict: AA - {}, NUM - {}".format(gamma_file_dict["AA"][i],gamma_file_dict["AA"][i][:-3]))
                gamma_file_dict["seq_num"].append(deepcopy(gamma_file_dict["AA"][i][:-3]))
    except:
        if gamma_file_dict["AA"][0][:-1].isnumeric():
            gamma_file_dict["seq_num"]=[gamma_file_dict["AA"][0][:-1]]
            for i in range(1,len(gamma_file_dict["AA"])):
                gamma_file_dict["seq_num"].append(deepcopy(gamma_file_dict["AA"][i][:-1]))

    return gamma_file_dict


def Read_peak_uncertainty(file_director:str, list_verson:str, s_dim:int, peak_objects:list[CResidue], version:int):
    Flag = False
    peaklistfile = f"{file_director}peak_lists/{list_verson}_peaks_noise.list"

    if os.path.exists(peaklistfile):
        RaportBox.write(f"\n\nLISTA:{peaklistfile}\n")
        with open(peaklistfile, 'r') as peak_noise_file:  
            p_lines_a = peak_noise_file.readlines()
            for indexl, line in enumerate(p_lines_a):
                if indexl > 0 and len(line) > 1:
                    item_a = line.split()
                    description = item_a[0]
                    aminoacids_number = read_aminoacid_number(description)
                    # aminoacids = description.split("-")
                    # aminoacids_number = ''.join([n for n in str(aminoacids[s_dim-1]) if n.isdigit()])
                    # if aminoacids_number == '':
                    #     aminoacids_number = ''.join([n for n in str(aminoacids[0]) if n.isdigit()])
                    RaportBox.write(f"aminoacids_number: {aminoacids_number,description}\n")
                    aminoacids_number = int(aminoacids_number)
                    for residue in peak_objects:
                        # residue = CResidue()
                        if aminoacids_number == residue.aa_number:
                            residue.peak_uncertainty[version] = float(item_a[-1])
                            RaportBox.write(f"peak_uncertainty = {float(item_a[-1])}\n")
        Flag = True
    else:
        Flag = False
        print_raport (f"There is no file like: {peaklistfile}!, so there will be used averange noise value\n")
    return Flag



"""Calc and write functions"""


def CheckOverlap(peak_list):
    # pl = peak_list
    for indexp1, peak1 in enumerate(peak_list):
        if peak1.is_peak:
            for indexp2, peak2 in enumerate(peak_list):
                if indexp2 < indexp1:
                    if peak2.is_peak:
                        peaks_distance = calc_distance(peak1.peak_pos_points, peak2.peak_pos_points)
                        if peaks_distance <= 3:
                            peak1.is_overlap = True
                            peak2.is_overlap = True
                            peak1.overlap_peaks[peak2.aa_number]=peaks_distance
                            peak2.overlap_peaks[peak1.aa_number]=peaks_distance
    return

def CalcErrorValue(peak_list,peak_number,angle_position,number_scan,Tc,uncertainty_val):
    Ix = peak_list[peak_number].peak_intens[1]
    Ia = peak_list[peak_number].peak_intens[0]
    NSa = number_scan[0]
    NSx = number_scan[1]
    if uncertainty_val == None:
        uncertainty_Ia = peak_list[peak_number].peak_uncertainty[0]
        uncertainty_Ix = peak_list[peak_number].peak_uncertainty[1]
        print_raport(f"{peak_list[peak_number].descript}: Peak_uncertainty for CalcErrorValue: {uncertainty_Ia}, {uncertainty_Ix}")
    else:
        uncertainty_Ia = uncertainty_val[0]
        uncertainty_Ix = uncertainty_val[1]
    squareDerivativeArcTanh = math.pow(1/(1-math.pow((Ix*NSa)/(Ia*NSx),2)),2)
    sumOfSquaresOfPeakHightDeviation = math.pow(uncertainty_Ix/Ia,2) + math.pow(uncertainty_Ia*Ix/math.pow(Ia,2),2)
    squareError = math.pow(1/Tc,2) * squareDerivativeArcTanh * math.pow(NSa/NSx,2) * sumOfSquaresOfPeakHightDeviation
    peak_list[peak_number+angle_position].ccrrate_error_value = math.sqrt(squareError)

    RaportBox.write("Calculating gamma error for peak number: {} \n squareDerivativeArcTanh = {} \n sumOfSquaresOfPeakHightDeviation = {}\n squareError = {}\n \t\t--> {}\n".format(peak_number,
                                                                                                                                                       squareDerivativeArcTanh,
                                                                                                                                                       sumOfSquaresOfPeakHightDeviation,
                                                                                                                                                       squareError,
                                                                                                                                                       math.sqrt(squareError)))



def CalcCCRRate(peak_list,angle_position,number_scan,Tc,CCR_name,noise_val,peak_uncertaintyflag):
    for indexp, peak in enumerate(peak_list):
        if not peak.is_overlap and peak.is_peak:
            peak_list[indexp+angle_position].is_ccr_rate = True
            if peak_uncertaintyflag:
                CalcErrorValue(peak_list,indexp,angle_position,number_scan,Tc,None)
            else:
                CalcErrorValue(peak_list,indexp,angle_position,number_scan,Tc,noise_val)
            try:
                # print (peak.peak_intens[1], number_scan[0], peak.peak_intens[0],number_scan[1], "---->", (peak.peak_intens[1]*number_scan[0])/(peak.peak_intens[0]*number_scan[1]))
                if CCR_name == "CCR_5" or CCR_name == "CCR_6":
                    ccr_rate_vol = -atanh((peak.peak_intens[1]*number_scan[0])/(peak.peak_intens[0]*number_scan[1]))/Tc
                else:
                    ccr_rate_vol = atanh((peak.peak_intens[1]*number_scan[0])/(peak.peak_intens[0]*number_scan[1]))/Tc
            except:
                ccr_rate_vol = (peak.peak_intens[1]*number_scan[0])/(peak.peak_intens[0]*number_scan[1]) #"atanh(x) - x shoudl be beetween -1 and 1"
                peak.ccrrate_calculation_error = True
            peak_list[indexp+angle_position].ccr_rate = ccr_rate_vol
            
    return

def Additional_text(exp):
    if exp.other=="":
        add_text=""
    else: add_text="_"+exp.other
    return add_text

def Add_ref_gamma(ref_seq_num,ref_CCR,peaks):
    # print ("ref_seq_num",ref_seq_num)
    for one_peak in peaks:
        for indexr, r_seq_num in enumerate(ref_seq_num):
            # print ("aa_number: {}, seq_number: {}".format(one_peak.aa_number,r_seq_num))
            if int(one_peak.aa_number) == int(r_seq_num) and ref_CCR[indexr]!='nan':
                one_peak.is_gamma_calc = True
                one_peak.gamma_ref = float(ref_CCR[indexr])
                # print ("gamma_ref",one_peak.gamma_ref)
    
    return

def calc_theor_Ix(peaks, angle_position, Tctime, number_scan, CCR_name):
    
    for indexp, one_peak in enumerate(peaks):
        
        if len(peaks)> indexp-angle_position >= 0 :
            residue = str(peaks[indexp+angle_position].aa_number)+peaks[indexp+angle_position].aa_name
            print_raport(residue)
            auto_peak_intensity = peaks[indexp-angle_position].peak_intens[0]
            if one_peak.is_gamma_calc and auto_peak_intensity != 0.0:
                peaks[indexp-angle_position].Ix_theor = float(math.tanh(one_peak.gamma_ref*Tctime)*auto_peak_intensity*number_scan[1]/number_scan[0])
                if CCR_name == "CCR_5" or CCR_name == "CCR_6":
                    peaks[indexp-angle_position].Ix_theor = -1 * peaks[indexp-angle_position].Ix_theor
                peaks[indexp-angle_position].Ix_theor_Ia_ratio = one_peak.Ix_theor/auto_peak_intensity
                peaks[indexp-angle_position].Ix_theor_Ia_ratio_without_NS = one_peak.Ix_theor_Ia_ratio*number_scan[0]/number_scan[1]
                
            else:
                residue = str(peaks[indexp+angle_position].aa_number)+peaks[indexp+angle_position].aa_name
                print_raport("This residue ({}) don't have reference value of gamma or referens (auto) peak".format(residue))


def LRegression_expresion(x,y):
    slope, intercept, r, p, std_err = stats.linregress(x, y)
    matching_factor_for_y_x = (1-slope)/std_err #if it is <1 is good, >1 means that we have weak corellation x=y 
    r2 = r**2
    lr_exp = '{:.3f}*x + {:.3f}'.format(slope,intercept)
    RaportBox.write(f"\nlinear regretion: {lr_exp}, r2: {r**2:.3f}\n p: {p}, std_err:{std_err}\n")
    lr_exp_y_list = []
    for i in x:
        lr_exp_y_list.append(deepcopy(slope*i+intercept))
    return lr_exp, r2, slope, intercept, matching_factor_for_y_x

def WeightedLRegression_expresion(x,y,uncertainty_val):
    regr = LinearRegression()
    # X = np.array([[i,i] for i in x if x!='nan'])
    # Y = np.array([[i,y[index]] for index, i in enumerate(x)])
    X = np.array([[i] for i in x if x!='nan'])
    Y = np.array([[i] for i in y])
    # X.reshape(-1,1)
    # Y.reshape(-1,1)
    
    weight_val = [1/i for i in uncertainty_val]
    regr.fit(X, Y, weight_val)
    predict_y = [i for [i] in regr.predict(X).tolist()]
    
    r2 = regr.score(X, Y, weight_val)
    [[slope]] = regr.coef_
    intercept = regr.intercept_[0] #type: ignore
    RaportBox.write(f"slope: {slope:.3f}, intercept: {intercept:.3f}")
    # print(f"r2: {r2:.3f}, \nregr.coef_:{regr.coef_},\nreg.intercept_:{regr.intercept_}")
    lr_exp = '{:.3f}*x + {:.3f}'.format(slope,intercept)
    RaportBox.write(f"weighted linear regretion: {lr_exp}, r2: {r2:.3f}")
    return lr_exp, r2,slope, intercept

def calc_R2(exp_list,theory_list):
    SSR = 0.0
    SSE = 0.0
    SST = 0.0
    average_exp = round(sum(exp_list)/len(exp_list), 3)
    for indext, theory_val in enumerate(theory_list):
        SSR+=(theory_val-average_exp)**2
        SSE+=(exp_list[indext]-theory_val)**2
        SST+=(exp_list[indext]-average_exp)**2
    R2_val = SSR/SST
    R2_val2 = 1 - (SSE/SST)
    print_raport ("SSR = {:.0f}; SSE = {:.0f}; SST = {:.0f}; SSR+SSE = {:.0f}; R2_1 = {:.4f}; R2_2 = {:.4f}".format(SSR, SSE, SST, SSR+SSE, R2_val, R2_val2))
    return R2_val2



def WeightedLRegression_expresion_by_hand(x,y,uncertainty_val):
    RaportBox.write("\nby Hand\n")
    list_of_invert_squer_uncertainty = []
    list_of_avg_x = []
    list_of_avg_square_x = []
    list_of_avg_y = []
    list_of_x_avgx_y_avgy_div_uncertainty = []
    list_of_x_avgx_squer_div_uncertainty = []


    for i in range(len(x)):
        list_of_invert_squer_uncertainty.append(deepcopy(1/(uncertainty_val[i]**2)))
        list_of_avg_x.append(deepcopy(x[i]/(uncertainty_val[i]**2)))
        list_of_avg_y.append(deepcopy(y[i]/(uncertainty_val[i]**2)))
        list_of_avg_square_x.append(deepcopy((x[i]**2)/(uncertainty_val[i]**2)))

    sum_of_invert_squer_uncertainty = sum(list_of_invert_squer_uncertainty)
    avg_x = sum(list_of_avg_x)/sum_of_invert_squer_uncertainty
    avg_y = sum(list_of_avg_y)/sum_of_invert_squer_uncertainty
    avg_square_x = sum(list_of_avg_square_x)/sum_of_invert_squer_uncertainty

    for i in range(len(x)):
        list_of_x_avgx_y_avgy_div_uncertainty.append(deepcopy((x[i]-avg_x)*(y[i]-avg_y)/(uncertainty_val[i]**2)))
        list_of_x_avgx_squer_div_uncertainty.append(deepcopy(((x[i]-avg_x)/uncertainty_val[i])**2))

    sum_of_x_avgx_y_avgy_div_uncertainty = sum(list_of_x_avgx_y_avgy_div_uncertainty)
    sum_of_x_avgx_squer_div_uncertainty = sum(list_of_x_avgx_squer_div_uncertainty)
    print_raport(f"sum_of_x_avgx_y_avgy_div_uncertainty: {sum_of_x_avgx_y_avgy_div_uncertainty}, sum_of_x_avgx_squer_div_uncertainty: {sum_of_x_avgx_squer_div_uncertainty}")
    slope = sum_of_x_avgx_y_avgy_div_uncertainty/sum_of_x_avgx_squer_div_uncertainty  
    intercept = avg_y-(slope*avg_x)
    slope_uncertainty = math.sqrt(1/(sum_of_x_avgx_squer_div_uncertainty))
    intercept_uncertainty = math.sqrt((avg_square_x)/sum_of_x_avgx_squer_div_uncertainty)

    print_raport(f"slope: {slope:.3f},intercept:{intercept:.3f},slope_uncertainty:{slope_uncertainty:.5f},intercept_uncertainty:{intercept_uncertainty:.5f}")
    lr_exp = '{:.2f}*x + {:.2f}'.format(slope,intercept)
    print_raport(f"weighted linear regretion: {lr_exp}")
    factor_a_for_y_x = (1-abs(slope))/slope_uncertainty
    factor_b_for_y_x = (1-abs(intercept))/intercept_uncertainty
    print_raport(f"1-a/delta_a = {factor_a_for_y_x}")
    y_from_line = []
    for one_x in x:
        y_from_line.append(deepcopy(one_x*slope+intercept))
    r2 = calc_R2(y,y_from_line)

    return lr_exp, slope, slope_uncertainty, intercept, intercept_uncertainty, factor_a_for_y_x, factor_b_for_y_x, r2


def check_if_min_max(suspect_min:Union[int,float], 
                     suspect_max:Union[int,float], 
                     min_max_list:list[Union[int,float]]) -> list[Union[int,float]]:
    if suspect_min<min_max_list[0]:
        min_max_list[0]=suspect_min
    if suspect_max<min_max_list[0]:
        min_max_list[0]=suspect_max
    if suspect_min>min_max_list[1]:
        min_max_list[1]=suspect_min
    if suspect_max>min_max_list[1]:
        min_max_list[1]=suspect_max
    return min_max_list
    
def setup_plot_area(plot, min_max_list:list[Union[int,float]]):
    additional_space = (abs(min_max_list[1])+abs(min_max_list[0]))/4
    plot.set_xlim(min_max_list[0]-additional_space,min_max_list[1]+additional_space)
    plot.set_ylim(min_max_list[0]-additional_space,min_max_list[1]+additional_space)
    plot.set_aspect('equal', adjustable='box')
    return plot

def plot_gamma_gamma(peaks,ccr_name,transparent_plot=False,add="",add2=""):

    min_max_value = [+100.0,-100.0]         #minimal and maximal value of CCR rate or reference gamma - it is necessary to make plot  
    gamma_calculated = []
    gamma_experimental = []
    gamma_calc_error = []            # ERROR GAMMA INFO, 9.11.2023
    gamma_differences = []
    for one_peak in peaks:
        if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
            # print(one_peak.gamma_ref,one_peak.ccr_rate)
            gamma_calculated.append(deepcopy(one_peak.gamma_ref))
            gamma_experimental.append(deepcopy(one_peak.ccr_rate))
            gamma_calc_error.append(deepcopy(one_peak.ccrrate_error_value))            # ERROR GAMMA INFO, 9.11.2023
            min_max_value = check_if_min_max(one_peak.gamma_ref,one_peak.ccr_rate,min_max_value)
            gamma_differences.append(deepcopy(one_peak.ccr_rate-one_peak.gamma_ref))

    plt.rcParams['font.size'] = '14'

    wilcoxon_rank = stats.wilcoxon(gamma_calculated,gamma_experimental)[1] # results: statistic, pvalue, zstatistic  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
    tStudent_2 = stats.ttest_ind(gamma_calculated,gamma_experimental)[0]    # results: statistic, pvalue, df  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
    # pvalue is for alternative hypothesis, so:
    # if p>a -> H0, 
    # if p<a -> H1
    tStudent_3 = stats.ttest_rel(gamma_calculated,gamma_experimental)[0] # results: statistic, pvalue, 
    avg_distance = np.mean(gamma_differences)
    std_distance = np.std(gamma_differences)
    t_critical_val = stats.t.ppf(q=1-.05,df=len(gamma_differences))
    tStudent_1 = avg_distance/(std_distance/math.sqrt(len(gamma_differences))) #type: ignore

    label = f'structure-predicted vs experimental,\nwilcoxon rank = {wilcoxon_rank:.4f}\ntStudent = {tStudent_2:.4f}\ntStudent for y-x = {tStudent_1:.4f}\nT critical val = {t_critical_val:.4f}\ntStudent_3 = {tStudent_3:.4f}'
    
    # Plotting both the curves simultaneously
    plt.axline([0,0],slope=1, linestyle=(0, (5, 5)), linewidth=1.5, color='darkgray', label='x=y')
    plt.scatter(gamma_calculated, gamma_experimental, s=5, color='#252525ff', 
                label=label)
    plt.errorbar(gamma_calculated, gamma_experimental, 
                    yerr=gamma_calc_error, 
                    fmt='none', color='#252525ff') #

    print_raport(f"\n {ccr_name}")
    print_raport (label)
    # print(f"gamma_calculated:{gamma_calculated} \ngamma_experimental: {gamma_experimental}\ngamma_calc_error: {gamma_calc_error}")
    linear_expretion,linear_r2, linear_slope, linear_intercept, matching_factor_for_y_x_LR= LRegression_expresion(gamma_calculated, gamma_experimental)
    weighted_expretion,weighted_r2, weighted_slope, weighted_intercept = WeightedLRegression_expresion(gamma_calculated, gamma_experimental,gamma_calc_error)
    weighted_expretion2, weighted_slope2, weighted_slope2error, weighted_intercept2, weighted_intercept2error, factor_a, factor_b, r2 = WeightedLRegression_expresion_by_hand(gamma_calculated, gamma_experimental,gamma_calc_error)
    plt.axline([0,linear_intercept],slope=linear_slope, linestyle=(0, (5, 5)), linewidth=1.5, color='red', label=f'l: {linear_expretion}, r2 = {linear_r2:.3f}, \nmatching factor for y=x: {matching_factor_for_y_x_LR:.1f}')
    plt.axline([0,weighted_intercept],slope=weighted_slope, linestyle=(0, (5, 5)), linewidth=1.5, color='blue', label=f'lw: {weighted_expretion}, r2 = {weighted_r2:.3f}')
    plt.axline([0,weighted_intercept2],slope=weighted_slope2, linestyle=(0, (5, 5)), linewidth=1.5, color='purple', label=f'lw2: {weighted_expretion2}, \nmatching factor for y=x: {factor_a:.1f}')

    plt.title(ccr_name+add)
    # Naming the x-axis, y-axis and the whole graph 
    plt.xlabel('structure-predicted \u0393, $s^{-1}$')
    plt.ylabel('experimental \u0393, $s^{-1}$')
    plt.xlim(min_max_value[0]-5.0,min_max_value[1]+5.0)
    plt.ylim(min_max_value[0]-5.0,min_max_value[1]+5.0)
    # plt.ticklabel_format(style='plain')
    start, end = plt.gca().get_ylim()
    start = (start//5)*5
    end = ((end//5)+1)*5
    if abs(start)+abs(end)<= 30:
        plt.gca().yaxis.set_ticks(np.arange(int(start), int(end), 5))
        plt.gca().xaxis.set_ticks(np.arange(int(start), int(end), 5))
    else:
        plt.gca().yaxis.set_ticks(np.arange(int(start), int(end), 10))
        plt.gca().xaxis.set_ticks(np.arange(int(start), int(end), 10))

    plt.gca().set_aspect('equal', adjustable='box')
    # Adding legend, which helps us recognize the curve according to it's color
    plt.legend(fontsize="small",bbox_to_anchor=(1.05, 1),
                         loc='upper left', borderaxespad=0.)
    
    # To load the display window
    # plt.savefig("{}/wykresy2/{}_exp_vs_calc.png".format(Dir, tfn[0].name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
    plt.savefig("{}{}_exp_vs_calc{}.png".format(file_director, ccr_name,add2), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
    # plt.show()
    plt.clf()
    plt.close()




def prepare_gamma_calc_vs_gamma_exp_table(peaks_tab:list[list[CResidue]],
                                          min_max_value:list) -> tuple[list[list[float]],list[Union[int,float]]]:
    set_of_data = []
    
    for indext in range(len(peaks_tab)):
        gamma_calculated = []
        gamma_experimental = []
        gamma_calc_error = []
        for one_peak in peaks_tab[indext]:
            if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
                gamma_calculated.append(deepcopy(one_peak.gamma_ref))
                gamma_experimental.append(deepcopy(one_peak.ccr_rate))
                gamma_calc_error.append(deepcopy(one_peak.ccrrate_error_value))
                min_max_value = check_if_min_max(one_peak.gamma_ref,one_peak.ccr_rate,min_max_value)
        set_of_data.append(deepcopy([gamma_calculated,gamma_experimental,gamma_calc_error]))
    return set_of_data,min_max_value


def prepare_gamma_calc_vs_gamma_exp_table_diff_min_max(peaks_tab:list[list[CResidue]],) -> tuple[list[list[float]],list[list[Union[int,float]]]]:
    set_of_data = []
    set_of_min_max_value = []
    for indext in range(len(peaks_tab)):
        gamma_calculated = []
        gamma_experimental = []
        gamma_calc_error = []
        min_max_value = [+100.0,-100.0] #minimal and maximal value of CCR rate or reference gamma - it is necessary to make plot  
        for one_peak in peaks_tab[indext]:
            if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
                gamma_calculated.append(deepcopy(one_peak.gamma_ref))
                gamma_experimental.append(deepcopy(one_peak.ccr_rate))
                gamma_calc_error.append(deepcopy(one_peak.ccrrate_error_value))
                min_max_value = check_if_min_max(one_peak.gamma_ref,one_peak.ccr_rate,min_max_value)
        set_of_data.append(deepcopy([gamma_calculated, gamma_experimental, gamma_calc_error]))
        set_of_min_max_value.append(deepcopy(min_max_value))
    return set_of_data,set_of_min_max_value


def prepare_gamma_exp_vs_diff_nus_table(peaks_tab:list[list[CResidue]],
                                        min_max_value:list,descrition:list[str]) -> tuple[list[list[float]],list[str],list[Union[int,float]]]:
    set_of_data = [[],[]]
    
    max_number = 0
    min_number = 10000000
    min_exp = -1
    max_exp = -1 
    for indext, text in enumerate(descrition):
        nus_number = text[1:-3]
        if nus_number.isdigit() == True:
            nus_number = int(nus_number)
            if nus_number > max_number:
                max_number = nus_number
                max_exp = indext
            if nus_number < min_number:
                min_number = nus_number
                min_exp = indext
    min_max_discript = [descrition[min_exp],descrition[max_exp]]
    peaks_table_min = peaks_tab[min_exp]
    peaks_table_max = peaks_tab[max_exp]
    for indexp, one_peak in enumerate(peaks_table_max):
        if one_peak.is_ccr_rate and peaks_table_min[indexp].is_ccr_rate and one_peak.aa_name!="G":
            set_of_data[0].append(deepcopy(peaks_table_min[indexp].ccr_rate))
            set_of_data[1].append(deepcopy(one_peak.ccr_rate))
            min_max_value = check_if_min_max(peaks_table_min[indexp].ccr_rate,one_peak.ccr_rate,min_max_value)
        
    return set_of_data,min_max_discript,min_max_value


def plot_gamma_gamma_all_together(peaks_tab:list[list[CResidue]],ccr_name:list[str],add:list[str],transparent_plot=False, nrows=3, ncols=3, add2="", publication_style=False):

    
    fig, axs = plt.subplots(nrows=nrows,ncols=ncols)
    if not publication_style:
        fig.suptitle(f"Comparision of experimental CCR rates with structure-predicted CCR rates based on:\n\"{gamma_cal_file_name}\" ")
    
    set_of_data_all_type_of_NUS,min_max_value = prepare_gamma_calc_vs_gamma_exp_table_diff_min_max(peaks_tab)
    
   # first row - plot each NUS number in particular plot
    for indexa, ax in enumerate(axs.flat):
        if indexa < len(peaks_tab):
            ax.axline([0,0],slope=1, linestyle=(0, (3, 3)), linewidth=1, color='darkgray') 
            ax.scatter(set_of_data_all_type_of_NUS[indexa][0],set_of_data_all_type_of_NUS[indexa][1],s=3, color='#252525ff') #
            ax.errorbar(set_of_data_all_type_of_NUS[indexa][0], set_of_data_all_type_of_NUS[indexa][1], 
                    yerr=set_of_data_all_type_of_NUS[indexa][2], 
                    fmt='none', color='#252525ff') #
            weighted_expretion2, weighted_slope2, weighted_slope2error, weighted_intercept2, weighted_intercept2error, factor_a, factor_b, r2 = WeightedLRegression_expresion_by_hand(set_of_data_all_type_of_NUS[indexa][0],
                                                                                                                                       set_of_data_all_type_of_NUS[indexa][1],
                                                                                                                                       set_of_data_all_type_of_NUS[indexa][2])
            ax.axline([0,weighted_intercept2],slope=weighted_slope2, linestyle=(0, (5, 5)), linewidth=2, 
                      color='mediumorchid', label=f'{add[indexa][1:]}:    {len(set_of_data_all_type_of_NUS[indexa][1])} peaks,    {weighted_expretion2},  delta_a: {weighted_slope2error:.3f}, delta_b: {weighted_intercept2error:.3f}    r2 = {r2:.2f}    factor \'a\' = {factor_a:.1f}, factor \'b\' = {factor_b:.1f}')
            ax.set_title(f"{ccr_name[indexa]} {add[indexa][1:]}", fontsize=10) #
            ax = setup_plot_area(ax,min_max_value[indexa])
            ax.tick_params(labelsize=8)

    #mediumpurple
    plt.subplots_adjust(hspace=0.3,wspace = 0.3)

    if add2 == '':
        add2 = add[0]

    if publication_style:
        fig.supxlabel('structure-predicted \u0393, $s^{-1}$', fontsize=15)
        fig.supylabel('experimental \u0393, $s^{-1}$', fontsize=15)
        cm = 1/2.54  # centimeters in inches
        fig.set_figwidth(17.4*cm)
        fig.set_figheight(17*cm)
        fig.savefig("{}{}_all_exp_vs_calc{}.png".format(file_director,gamma_cal_file_name[2:-4],add2), bbox_inches="tight", 
                pad_inches=0.3, transparent=transparent_plot)
        fig.savefig("{}{}_all_exp_vs_calc{}.svg".format(file_director,gamma_cal_file_name[2:-4],add2), bbox_inches="tight", 
                pad_inches=0.3, transparent=transparent_plot)
    else:
        fig.supxlabel('structure-predicted \u0393, $s^{-1}$')
        fig.supylabel('experimental \u0393, $s^{-1}$')
        fig.set_figwidth(ncols*5)
        fig.set_figheight(nrows*5)
        fig.savefig("{}_all_exp_vs_calc{}.png".format(file_director,add2), bbox_inches="tight", 
                pad_inches=0.3, transparent=transparent_plot)
    plt.close()
    plt.clf()




def plot_gamma_gamma_and_diff_together(peaks_tab:list[list[CResidue]],ccr_name:str,add:list[str],transparent_plot=False,publication_style=False):

    # plt.rcParams['font.size'] = '14'

    min_max_value = [+100.0,-100.0]
    fig, axs = plt.subplots(nrows=2,ncols=len(peaks_tab))
    if not publication_style:
        fig.suptitle(f"Comparision of CCR rates ({ccr_name}) for diffrent experiment variant with structure-predicted CCR rates based on:\n\"{gamma_cal_file_name}\"")
    
    set_of_data_all_type_of_NUS,min_max_value = prepare_gamma_calc_vs_gamma_exp_table(peaks_tab,min_max_value)

    set_of_data_min_max_NUS,min_max_discript,min_max_value = prepare_gamma_exp_vs_diff_nus_table(peaks_tab,min_max_value,add)
   
   # first row - plot each NUS number in particular plot
    for indext in range(len(peaks_tab)):
        axs[0,indext].axline([0,0],slope=1, linestyle=(0, (3, 3)), linewidth=1, color='darkgray') 
        axs[0,indext].scatter(set_of_data_all_type_of_NUS[indext][0],set_of_data_all_type_of_NUS[indext][1],s=10,) #
        axs[0,indext].errorbar(set_of_data_all_type_of_NUS[indext][0], set_of_data_all_type_of_NUS[indext][1], 
                    yerr=set_of_data_all_type_of_NUS[indext][2], 
                    fmt='none', color='#252525ff') #
        if len(set_of_data_all_type_of_NUS[indext][0])>0:
            weighted_expretion2, weighted_slope2, weighted_slope2error, weighted_intercept2, weighted_intercept2error, factor_a, factor_b, r2 = WeightedLRegression_expresion_by_hand(set_of_data_all_type_of_NUS[indext][0],
                                                                                                                                       set_of_data_all_type_of_NUS[indext][1],
                                                                                                                                       set_of_data_all_type_of_NUS[indext][2])
            axs[0,indext].axline([0,weighted_intercept2],slope=weighted_slope2, linestyle=(0, (5, 5)), linewidth=1.5, 
                      color='purple', label=f'{add[indext][1:]}:    {len(set_of_data_all_type_of_NUS[indext][1])} peaks,    {weighted_expretion2},  delta_a: {weighted_slope2error:.3f}, delta_b: {weighted_intercept2error:.3f}\nr2 = {r2:.2f}    factor \'a\' = {factor_a:.1f}, factor \'b\' = {factor_b:.1f}')
            # axs[0,indext].legend()
        axs[0,indext].set_title(add[indext][1:], fontsize=10)
        axs[0,indext] = setup_plot_area(axs[0,indext],min_max_value)
        axs[0,indext].tick_params(labelsize=8)


    # second row
    axs[1,0].axline([0,0],slope=1, linestyle=(0, (3, 3)), linewidth=1, color='darkgray', label='x=y') 
    for indext in range(len(peaks_tab)):
        axs[1,0].scatter(set_of_data_all_type_of_NUS[indext][0],set_of_data_all_type_of_NUS[indext][1],
                         s=10,label=str(add[indext][1:])) #
    axs[1,0] = setup_plot_area(axs[1,0],min_max_value)
    axs[1,0].tick_params(labelsize=8)
    fig.legend(fontsize="medium",loc='lower left', bbox_to_anchor=(0.3, 0.1))

    for indext in range(1,len(peaks_tab)-1):
        axs[1][indext].set_axis_off()

    axs[1,len(peaks_tab)-1].axline([0,0],slope=1, linestyle=(0, (3, 3)), linewidth=1, color='darkgray') 
    axs[1,len(peaks_tab)-1].scatter(set_of_data_min_max_NUS[1],set_of_data_min_max_NUS[0],s=10,) #
    axs[1,len(peaks_tab)-1]= setup_plot_area(axs[1,len(peaks_tab)-1],min_max_value)
    axs[1,len(peaks_tab)-1].set_xlabel(min_max_discript[1]+' NUS \u0393, $s^{-1}$',fontsize=10)
    axs[1,len(peaks_tab)-1].set_ylabel(min_max_discript[0]+' NUS \u0393, $s^{-1}$',fontsize=10)
    axs[1,len(peaks_tab)-1].tick_params(labelsize=8)

    fig.supxlabel('structure-predicted \u0393, $s^{-1}$')
    fig.supylabel('experimental \u0393, $s^{-1}$')

    plt.subplots_adjust(hspace=0.3,wspace = 0.3)

    if publication_style:
        cm = 1/2.54  # centimeters in inches
        fig.set_figwidth(17.4*cm)
        fig.set_figheight(17*cm)
        fig.savefig("{}{}_exp_vs_calc_diff.svg".format(file_director, ccr_name), bbox_inches="tight", 
                pad_inches=0.3, transparent=transparent_plot)
    else:
        fig.set_figwidth(len(peaks_tab)*3.5)
        fig.set_figheight(7)
        fig.savefig("{}{}_exp_vs_calc_diff.png".format(file_director, ccr_name), bbox_inches="tight", 
                pad_inches=0.3, transparent=transparent_plot)
    plt.close()
    plt.clf()






def plot_error_histogram(peaks,ccr_name,transparent_plot=False,add="",add2=""): # ERROR GAMMA INFO, 9.11.2023

    gamma_experimental = []
    gamma_calc_error = []            
    peak_intens_auto = []         
    peak_intens_cross = []
    for one_peak in peaks:
        if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
            gamma_experimental.append(deepcopy(one_peak.ccr_rate))
            gamma_calc_error.append(deepcopy(one_peak.ccrrate_error_value)) 
            peak_intens_auto.append(deepcopy(one_peak.peak_intens[0]))
            peak_intens_cross.append(deepcopy(one_peak.peak_intens[1]))

    # print ("gamma_calculated:\t", len(gamma_calculated))
    # print ("gamma_experimental:\t", len(gamma_experimental))
    # print ("gamma_calc_error:\t", gamma_calc_error)
    # print ("min_max_value:\t", min_max_value)
    fig, axs = plt.subplots(nrows=3)
    plt.subplots_adjust(hspace=0.4)
    fig.set_figheight(14)

    axs[0].scatter(gamma_experimental, gamma_calc_error, s=5, color='#0066ffff')
    axs[0].set_title("{} {}\n Experimental \u0393 vs Uncertainty value".format(ccr_name,add))
    axs[0].set_ylabel('Uncertainty value')
    axs[0].set_xlabel('Experimental \u0393, $s^{-1}$')

    axs[1].scatter(peak_intens_auto, gamma_calc_error, s=5, color='#0066ffff')
    axs[1].set_title("Peak intens in auto vs Uncertainty value")
    axs[1].set_ylabel('Uncertainty value')
    axs[1].set_xlabel('Peak intens in auto')

    axs[2].scatter(peak_intens_cross, gamma_calc_error, s=5, color='#0066ffff')
    axs[2].set_title("Peak intens in cross vs Uncertainty value")
    axs[2].set_ylabel('Uncertainty value')
    axs[2].set_xlabel('Peak intens in cross')

    # plt.title(ccr_name+add)
    plt.savefig("{}Hist_{}_exp_vs_Uncertainty{}.png".format(file_director, ccr_name,add2), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
    # plt.show()
    plt.clf()
    plt.close()





def plot_gamma_diff(peaks_tab,ccr_name,add,transparent_plot=False):

    plt.rcParams['font.size'] = '14'
    # Plotting both the curves simultaneously
    plt.axline([0,0],slope=1, linestyle=(0, (5, 5)), linewidth=1.5, color='darkgray', label='x=y')

    min_max_value = [+100.0,-100.0]         #minimal and maximal value of CCR rate or reference gamma - it is necessary to make plot  
    for indexp, peaks in enumerate(peaks_tab):
        gamma_calculated = []
        gamma_experimental = []
        for one_peak in peaks:
            if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
                gamma_calculated.append(deepcopy(one_peak.gamma_ref))
                gamma_experimental.append(deepcopy(one_peak.ccr_rate))
                if one_peak.gamma_ref-one_peak.ccr_rate>5 and indexp==0:
                    plt.annotate(one_peak.aa_number,(one_peak.gamma_ref,one_peak.ccr_rate),
                                 textcoords="offset points",xytext=(0,-12),ha='center')
                if one_peak.ccr_rate-one_peak.gamma_ref>5 and indexp==0:
                    plt.annotate(one_peak.aa_number,(one_peak.gamma_ref,one_peak.ccr_rate),
                                 textcoords="offset points",xytext=(0,10),ha='center')
                min_max_value = check_if_min_max(one_peak.gamma_ref,one_peak.ccr_rate,min_max_value)
        plt.scatter(gamma_calculated, gamma_experimental, s=10, label=str(add[indexp][1:]))

    plt.title("comparision of CCR rates for diffrent number of NUS points", fontsize=10)
    plt.suptitle(ccr_name)
    # Naming the x-axis, y-axis and the whole graph 
    plt.xlabel('structure-predicted \u0393, $s^{-1}$')
    plt.ylabel('experimental \u0393, $s^{-1}$')
    plt.xlim(min_max_value[0]-5.0,min_max_value[1]+5.0)
    plt.ylim(min_max_value[0]-5.0,min_max_value[1]+5.0)
    # plt.ticklabel_format(style='plain')
    start, end = plt.gca().get_ylim()
    start = (start//5)*5
    end = ((end//5)+1)*5
    if abs(start)+abs(end)<= 30:
        plt.gca().yaxis.set_ticks(np.arange(int(start), int(end), 5))
        plt.gca().xaxis.set_ticks(np.arange(int(start), int(end), 5))
    else:
        plt.gca().yaxis.set_ticks(np.arange(int(start), int(end), 10))
        plt.gca().xaxis.set_ticks(np.arange(int(start), int(end), 10))

    plt.gca().set_aspect('equal', adjustable='box')
    # Adding legend, which helps us recognize the curve according to it's color
    plt.legend(fontsize="small")
    
    # To load the display window
    # plt.savefig("{}/wykresy2/{}_exp_vs_calc.png".format(Dir, tfn[0].name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
    plt.savefig("{}{}_exp_vs_calc_diff.png".format(file_director, ccr_name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
    # plt.show()
    plt.clf()
    plt.close()





def plot_cross_intens_theor_exp_together(peaks_tab,ccr_name,add,transparent_plot=False):

    # plt.rcParams['font.size'] = '14'

    min_max_value = [+100.0,-100.0]
    fig, axs = plt.subplots(ncols=len((peaks_tab)))
    fig.suptitle("Comparision of intensity of peaks in transfer spectra {}\n(structure-predicted vs experimental) \nwith diffrent number of NUS points".format(ccr_name))
    
    set_of_data = []
    for indext in range(len(peaks_tab)):
        Ix_theory = []
        Ix_experimental = []
        for one_peak in peaks_tab[indext]:
            if one_peak.is_peak and one_peak.Ix_theor!=0.0:
                # print ("we have peak {} for spectra {} ({})".format(one_peak.aa_number,ccr_name,add[indexa]))
                # print("Intensity of Ix_theory: {}, Ix_experimental: {}, min_max_value: {}".format(one_peak.Ix_theor, one_peak.peak_intens[1],min_max_value))
                Ix_theory.append(deepcopy(one_peak.Ix_theor))
                Ix_experimental.append(deepcopy(one_peak.peak_intens[1]))
                min_max_value = check_if_min_max(one_peak.Ix_theor,one_peak.peak_intens[1],min_max_value)
        set_of_data.append(deepcopy([Ix_theory,Ix_experimental]))
        # print("lenth of Ix_theory: {}, Ix_experimental: {} for spectra {} ({})".format(len(Ix_theory), len(Ix_experimental),ccr_name,add[indext]))
        # print(Ix_theory, Ix_experimental)

    for indext in range(len(peaks_tab)):
        axs[indext].axline([0,0],slope=1, linestyle=(0, (3, 3)), linewidth=1, color='darkgray', label='x=y') 
        axs[indext].scatter(set_of_data[indext][0],set_of_data[indext][1],s=10,) #
        axs[indext].set_title(add[indext][1:],fontsize=14)
        axs[indext] = setup_plot_area(axs[indext],min_max_value)
        axs[indext].tick_params(labelsize=8)
        axs[indext].set_aspect('equal', adjustable='box')
    fig.supxlabel('structure-predicted intensity')
    fig.supylabel('experimental intensity')
    fig.set_figwidth(len(peaks_tab)*3)
    fig.set_figheight(4)
    # plt.subplots_adjust(hspace=0.3,wspace = 0.3)
    fig.savefig("{}{}_exp_vs_calc_Ix.png".format(file_director, ccr_name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
    plt.close()
    plt.clf()



def plot_error_histogram_together(peaks_tab,ccr_name,add,transparent_plot=False): # ERROR GAMMA INFO, 9.11.2023

    set_of_data = []
    for indext in range(len(peaks_tab)):
        gamma_experimental = []
        gamma_calc_error = []            
        peak_intens_auto = []         
        peak_intens_cross = []
        for one_peak in peaks_tab[indext]:
            if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
                gamma_experimental.append(deepcopy(one_peak.ccr_rate))
                gamma_calc_error.append(deepcopy(one_peak.ccrrate_error_value)) 
                peak_intens_auto.append(deepcopy(one_peak.peak_intens[0]))
                peak_intens_cross.append(deepcopy(one_peak.peak_intens[1]))
        set_of_data.append(deepcopy({"gamma_experimental":gamma_experimental,
                                      "gamma_calc_error":gamma_calc_error,
                                      "peak_intens_auto":peak_intens_auto,
                                      "peak_intens_cross":peak_intens_cross}))

    fig, axs = plt.subplots(nrows=3,ncols=len((peaks_tab)))    
    plt.subplots_adjust(hspace=0.4)
    fig.set_figheight(15)
    fig.set_figwidth(6*len((peaks_tab)))
    for indext in range(len(peaks_tab)):
        axs[0,indext].scatter(set_of_data[indext]["gamma_experimental"], set_of_data[indext]["gamma_calc_error"], s=5, color='#0066ffff')
        axs[0,indext].set_title("{} {}\n Experimental \u0393 vs Uncertainty value".format(ccr_name,str(add[indext][1:])))
        axs[0,indext].set_ylabel('Uncertainty value')
        axs[0,indext].set_xlabel('Experimental \u0393, $s^{-1}$')

        axs[1,indext].scatter(set_of_data[indext]["peak_intens_auto"], set_of_data[indext]["gamma_calc_error"], s=5, color='#0066ffff')
        axs[1,indext].set_title("Peak intens in auto vs Uncertainty value")
        axs[1,indext].set_ylabel('Uncertainty value')
        axs[1,indext].set_xlabel('Peak intens in auto')

        axs[2,indext].scatter(set_of_data[indext]["peak_intens_cross"], set_of_data[indext]["gamma_calc_error"], s=5, color='#0066ffff')
        axs[2,indext].set_title("Peak intens in cross vs Uncertainty value")
        axs[2,indext].set_ylabel('Uncertainty value')
        axs[2,indext].set_xlabel('Peak intens in cross')

    plt.savefig("{}Hist_{}_exp_vs_Uncertainty_all.png".format(file_director, ccr_name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
    # plt.show()
    plt.clf()
    plt.close()





def plot_gamma_exp_vs_epx(peaks_min,peaks_max,ccr_name,min_add,max_add,transparent_plot=False):

    plt.rcParams['font.size'] = '14'
    # Plotting both the curves simultaneously
    plt.axline([0,0],slope=1, linestyle=(0, (5, 5)), linewidth=1.5, color='darkgray', label='x=y')

    min_max_value = [+100.0,-100.0]         #minimal and maximal value of CCR rate or reference gamma - it is necessary to make plot  
    
    gamma_min = []
    gamma_max = []
    for indexp, one_peak in enumerate(peaks_max):
        if one_peak.is_ccr_rate and peaks_min[indexp].is_ccr_rate and one_peak.aa_name!="G":
            gamma_min.append(deepcopy(peaks_min[indexp].ccr_rate))
            gamma_max.append(deepcopy(one_peak.ccr_rate))
            if peaks_min[indexp].ccr_rate-one_peak.ccr_rate>2:
                plt.annotate(one_peak.aa_number,(one_peak.ccr_rate,peaks_min[indexp].ccr_rate),textcoords="offset points",xytext=(0,10),ha='center')
            if one_peak.ccr_rate-peaks_min[indexp].ccr_rate>2:
                plt.annotate(one_peak.aa_number,(one_peak.ccr_rate,peaks_min[indexp].ccr_rate),textcoords="offset points",xytext=(0,-15),ha='center')
            min_max_value = check_if_min_max(peaks_min[indexp].ccr_rate,one_peak.ccr_rate,min_max_value)
        plt.scatter(gamma_max,gamma_min, s=10)

    plt.title("comparision of CCR rates for diffrent number of NUS points", fontsize=10)
    plt.suptitle(ccr_name)
    # Naming the x-axis, y-axis and the whole graph 
    plt.xlabel(max_add+' NUS \u0393, $s^{-1}$')
    plt.ylabel(min_add+' NUS \u0393, $s^{-1}$')
    # plt.xlabel('min NUS \u0393, $s^{-1}$')
    # plt.ylabel('max NUS \u0393, $s^{-1}$')
    
    plt.xlim(min_max_value[0]-5.0,min_max_value[1]+5.0)
    plt.ylim(min_max_value[0]-5.0,min_max_value[1]+5.0)
    # plt.ticklabel_format(style='plain')
    start, end = plt.gca().get_ylim()
    start = (start//5)*5
    end = ((end//5)+1)*5
    if abs(start)+abs(end)<= 30:
        plt.gca().yaxis.set_ticks(np.arange(int(start), int(end), 5))
        plt.gca().xaxis.set_ticks(np.arange(int(start), int(end), 5))
    else:
        plt.gca().yaxis.set_ticks(np.arange(int(start), int(end), 10))
        plt.gca().xaxis.set_ticks(np.arange(int(start), int(end), 10))

    plt.gca().set_aspect('equal', adjustable='box')
    # Adding legend, which helps us recognize the curve according to it's color
    # plt.legend(fontsize="small")
    
    # To load the display window
    # plt.savefig("{}/wykresy2/{}_exp_vs_calc.png".format(Dir, tfn[0].name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
    # plt.show()
    plt.savefig("{}{}_minNUS_vs_maxNUS.png".format(file_director, ccr_name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
    plt.clf()
    plt.close()





def WriteCCRRate_small(peak_list, ccr_name, add=""):
    new_list = "{}{}{}.csv".format(RaportDir,CCRname2Ratename(ccr_name),add)
    with open(new_list, mode='w', newline='') as csv_file:
        headers = ['AA','CCR rate']
        writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
        for one_peak in peak_list:
            one_row = {}
            one_row["AA"] = one_peak.aa_number
            if one_peak.is_ccr_rate == True and one_peak.ccrrate_calculation_error == False and one_peak.aa_name != "G":
                one_row["CCR rate"] = '{:.4f}'.format(one_peak.ccr_rate)
            else:
                one_row["CCR rate"] = "nan"
            writer.writerow(one_row)
    # print ("CCR rate for CCR efect between {} rates are in file: {}".format(CCRname2Ratename(ccr_name), new_list), file=RaportBox)
    RaportBox.write("CCR rate for CCR efect between {} rates are in file: {}\n".format(CCRname2Ratename(ccr_name), new_list))
    return




def WriteCCRRate_all_info(peak_list, ccr_name, s_dim, autoname, crossname, Tctime, number_scan,add=""):
    new_list = "{}{}{}_CCRrate.list".format(RaportDir,ccr_name,add)
    # print ("new_peak_list", new_list, file=RaportBox)
    max_lenth_discrip = 0
    max_lenth_intens = 0

    for p in peak_list:
        pd = str(p.aa_number)+p.aa_name
        if len(pd)>max_lenth_discrip:
            max_lenth_discrip=len(pd)
        if len(str(p.peak_intens[0]))>max_lenth_intens:
            max_lenth_intens=len(str(p.peak_intens[0]))
    
    with open(new_list, 'w') as listfile:
        print ("\tauto list name = {}\tcross list name = {}\tTc = {}".format(autoname, crossname, Tctime), file=listfile) 
        print ("{}\t{:{sentence_len1}}\t{:{sentence_len2}}\t{:{sentence_len2}}\t CCR rate \t CCR Uncertainty \t CCR theor \t Ix theor \t Ix/Ia".format("AA","peak position","intensity in auto","intensity in cross",sentence_len1=7*s_dim, sentence_len2=max_lenth_intens), file=listfile)
        for one_peak in peak_list:
            dict_line={}
            peak_descr = str(one_peak.aa_number)+one_peak.aa_name
            dict_line['AA'] = "\n{:{sentence_len}}".format(peak_descr, sentence_len=max_lenth_discrip)
            # print ("{:{sentence_len}}".format(peak_descr, sentence_len=max_lenth_discrip,), end="\t", file=listfile)
            if one_peak.is_peak == True:
                peak_pos = ""
                for i in range(s_dim):
                    peak_pos += "{:.3f}\t".format(one_peak.peak_pos[i])
                dict_line['peak_pos'] = peak_pos
                dict_line['Ia'] = "{:{sentence_len}}".format("{:.2e}".format(one_peak.peak_intens[0]), sentence_len=max_lenth_intens) 
                dict_line['Ix'] = "{:{sentence_len}}".format("{:.2e}".format(one_peak.peak_intens[1]), sentence_len=max_lenth_intens) 
            else:
                if one_peak.is_ccr_rate == True:
                    dict_line['peak_pos'] = "{:{sentence_len}}\t\t".format("Info from other peak", sentence_len=7*s_dim+max_lenth_intens*2)
                else:
                    dict_line['peak_pos'] = "{:{sentence_len}}\t\t".format("No visible peak", sentence_len=7*s_dim+max_lenth_intens*2) 
            

            if one_peak.aa_name == "G":
                dict_line['ccr_rate'] = "{:{sentence_len}}".format("glycine", sentence_len=18)   #21 digit, glycine - no CCR rate
            elif one_peak.is_overlap == True:
                dict_line['ccr_rate'] = "This peak overlap with: {}".format(one_peak.overlap_peaks)
            elif one_peak.ccrrate_calculation_error == False and one_peak.is_ccr_rate:
                dict_line['ccr_rate'] = '{:18.4f}'.format(one_peak.ccr_rate, sentence_len=18)
            elif one_peak.ccrrate_calculation_error == True:
                dict_line['ccr_rate'] = "peak intens error: Ix/Ia = {:8.4f}".format(one_peak.is_ccr_rate, sentence_len=18) 
            elif one_peak.is_peak == True and one_peak.is_ccr_rate == False:
                dict_line['ccr_rate'] = "{:{sentence_len}}".format("Info to other peak", sentence_len=18)
            else:
                dict_line['ccr_rate'] = " "*20
            

            if one_peak.ccrrate_error_value != 0.0:
                dict_line['Uncertainty'] = '{:8.4f}'.format(one_peak.ccrrate_error_value, sentence_len=10)            # ERROR GAMMA INFO, 9.11.2023
            else:
                dict_line['Uncertainty'] = " "*10  
            

            if one_peak.is_gamma_calc and one_peak.is_ccr_rate:
                dict_line['ccr_rate_theor'] = "{:8.4f}".format(one_peak.gamma_ref, sentence_len=10)
                dict_line['Ix_theor'] = "{:8.2e}".format(one_peak.Ix_theor, sentence_len=10)
                dict_line['Ix_Ia'] = "{:8.2f}".format(one_peak.Ix_theor_Ia_ratio, sentence_len=10)
            else:
                dict_line['ccr_rate_theor'] =  dict_line['Ix_theor'] = dict_line['Ix_Ia'] = " "*10  

            for item in dict_line:
                print (dict_line[item], end="\t", file=listfile)
    # print ("All info from peaks for {} rates are in file: {}".format(ccr_name, new_list), file=RaportBox)
    RaportBox.write("All info from peaks for {} rates are in file: {}\n".format(ccr_name, new_list))
    return


def WriteCCRRate_all_info_CSV(peak_list, ccr_name, s_dim, autoname, crossname, Tctime, add=""):
    new_list = "{}{}{}_CCRrate.csv".format(RaportDir,ccr_name,add)
    max_lenth_discrip = 0
    max_lenth_intens = 0
    for p in peak_list:
        pd = str(p.aa_number)+p.aa_name
        if len(pd)>max_lenth_discrip:
            max_lenth_discrip=len(pd)
        if len(str(p.peak_intens[0]))>max_lenth_intens:
            max_lenth_intens=len(str(p.peak_intens[0]))
    with open(new_list, mode='w', newline='') as csv_file:
        # headers = ['AA','peak position','intensity in auto','intensity in cross','CCR rate','Comments']
        headers = ['AA']
        for i in range(1,s_dim+1):
            headers.append(deepcopy('w{}'.format(i)))
        headers.extend(deepcopy(['intensity in auto','intensity in cross','CCR rate','Uncertainty','Comments','Reference Gamma', 'Ix Theor', 'Ix/Ia', 'Ix/Ia without NS']))
        writer_row = csv.writer(csv_file,delimiter=",")
        writer_row.writerow(['auto list name =',autoname])
        writer_row.writerow(['cross list name =',crossname])
        writer_row.writerow(['Tc =', Tctime])
        writer_row.writerow(['---','---','---','---','---','---','---','---','---'])
        # writer_row.writerow(['AA','{:^{sentence_len1}}'.format('peak position',sentence_len1=7*s_dim),'intensity in auto','intensity in cross','CCR rate','Comments'])
        writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
        writer.writeheader()
        for one_peak in peak_list:
            one_row = {}
            peak_descr = str(one_peak.aa_number)+one_peak.aa_name
            one_row["AA"] = peak_descr
            if one_peak.is_peak == True and one_peak.is_ccr_rate == True:
                # peak_position = str(one_peak.peak_pos[0])
                for i in range(s_dim):
                    one_row['w{}'.format(i+1)] = '{:.3f}'.format(one_peak.peak_pos[i])
                # one_row["peak position"] = peak_position
                if one_peak.aa_name == "G":
                    one_row["Comments"] = "glycine - no CCR rate"
                else:
                    one_row["intensity in auto"] = one_peak.peak_intens[0]
                    one_row["intensity in cross"] = one_peak.peak_intens[1]
                    if one_peak.ccrrate_calculation_error == False:
                        one_row["CCR rate"] = '{:.4f}'.format(one_peak.ccr_rate)
                    else:
                        one_row["Comments"] = "peak intens error"
                
            elif one_peak.is_peak == False and one_peak.is_ccr_rate == True: 
                if one_peak.aa_name == "G":
                    one_row["Comments"] = "glycine - no CCR rate"
                else:
                    one_row["CCR rate"] = '{:.4f}'.format(one_peak.ccr_rate)
                    one_row["Comments"] = "Info from other peak"
            elif one_peak.is_overlap == True:
                one_row["Comments"] = "This peak overlap with: {}".format(one_peak.overlap_peaks)
            elif one_peak.is_peak == True and one_peak.is_ccr_rate == False:
                for i in range(s_dim):
                    one_row['w{}'.format(i+1)] = '{:.3f}'.format(one_peak.peak_pos[i])
            else:
                one_row["Comments"] = "No visible peak"
            
            if one_peak.ccrrate_error_value != 0.0:
                one_row['Uncertainty'] = '{:.4f}'.format(one_peak.ccrrate_error_value)            # ERROR GAMMA INFO, 9.11.2023
            

            if one_peak.is_gamma_calc and one_peak.is_ccr_rate:
                one_row['Reference Gamma'] = "{:.4f}".format(one_peak.gamma_ref)
                one_row['Ix Theor'] = "{:.2e}".format(one_peak.Ix_theor)
                one_row['Ix/Ia'] = "{:.2f}".format(one_peak.Ix_theor_Ia_ratio)
                one_row['Ix/Ia without NS'] = "{:.2f}".format(one_peak.Ix_theor_Ia_ratio_without_NS)
            writer.writerow(one_row)
    # print ("All info from peaks for {} rates are in file: {}".format(ccr_name, new_list), file=RaportBox)
    RaportBox.write("All info from peaks for {} rates are in file: {}\n".format(ccr_name, new_list))
    return

def Write_ALL_CCRRate_CSV():
    new_list = "{0}CCRrate.csv".format(file_director)
    with open(new_list, mode='w', newline='') as csv_file:
        headers = ['AA']
        for experiment in Experiments:
            Add_textt=Additional_text(experiment)
            headers.append(deepcopy('{}{}'.format(experiment.CCR_name,Add_textt)))
            headers.append(deepcopy('{}_error'.format(experiment.CCR_name)))            # ERROR GAMMA INFO, 9.11.2023
        writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
        writer.writeheader()
        zero_row = {}
        zero_row['AA'] = 'res'
        for experiment in Experiments:
            Add_textt=Additional_text(experiment)
            zero_row['{}{}'.format(experiment.CCR_name,Add_textt)] = '{}'.format(CCRname2Ratename(experiment.CCR_name))
        writer.writerow(zero_row)
        for indexaa, aanumber in enumerate(Sequence_List):
            one_row = {}
            for experiment in Experiments:
                Add_textt=Additional_text(experiment)
                one_peak = experiment.peak[indexaa]
                one_row["AA"] = str(one_peak.aa_number)+Res1to3(one_peak.aa_name)
                if one_peak.is_ccr_rate == True and one_peak.ccrrate_calculation_error == False and one_peak.aa_name != "G":
                    one_row['{}{}'.format(experiment.CCR_name,Add_textt)] = '{:.4f}'.format(one_peak.ccr_rate)
                    one_row['{}_error'.format(experiment.CCR_name)] = '{:.4f}'.format(one_peak.ccrrate_error_value)            # ERROR GAMMA INFO, 9.11.2023

                else:
                    one_row['{}{}'.format(experiment.CCR_name,Add_textt)] = "nan"
                    one_row['{}_error'.format(experiment.CCR_name)] = "nan"            # ERROR GAMMA INFO, 9.11.2023
            writer.writerow(one_row)
    print ("All CCR rates are in:", new_list)
    return


def Write_zeroCCRrates_Analisis_CSV():
    new_list = "zeroCCRrates_Analisis.csv".format(RaportDir)
    
    s_dim = Experiments[0].n_dim
    error_residue_dict={}
    error_experiment_dict={}
    error_aa_dict = {}
    invisible_cross_peaks = {}
    aa_dict_nun = {}

    with open(new_list, mode='w', newline='') as csv_file:
        headers = ['CCR name','AA']
        for i in range(1,s_dim+1):
            headers.append(deepcopy('w{}'.format(i)))
        headers.extend(deepcopy(['intensity in auto','Noise in auto','intensity in cross','Noise in cross','CCR rate','Uncertainty','Comments','Reference Gamma', 'Ix Theor', 'Ix/Ia', 'Ix/Ia without NS']))
        writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
        writer.writeheader()
        for experiment in Experiments:
            CCR_name =  experiment.CCR_name+Additional_text(experiment)
            noise_level = experiment.noise
            if experiment.is_peaklist == True:
                writer.writerow({})
                invisible_cross_peaks[CCR_name]=[0,0]
                error_experiment_dict[CCR_name]=0
                for peak_number, one_peak in enumerate(experiment.peak):
                    if one_peak.aa_name in aa_dict_nun:
                        aa_dict_nun[one_peak.aa_name]+=1
                    else:
                        aa_dict_nun[one_peak.aa_name]=1                    
                    deepVisibilityAnalisis(one_peak,experiment,invisible_cross_peaks)

                    if one_peak.is_ccr_rate and one_peak.aa_name != "G" and one_peak.is_gamma_calc:
                        if (one_peak.gamma_ref - one_peak.ccr_rate)/one_peak.ccr_rate >=2 or (one_peak.gamma_ref - one_peak.ccr_rate)/one_peak.ccr_rate <= -2:
                            
                            row_list = deepResidueAnalisys(writer,one_peak,peak_number,experiment,
                                                    error_residue_dict,
                                                    error_experiment_dict,
                                                    error_aa_dict)

        writer_row = csv.writer(csv_file,delimiter=",")
        writer_row.writerow("")
        
        deepExpAnalisys(writer_row,
                 error_residue_dict,
                 error_experiment_dict,
                 error_aa_dict,
                 invisible_cross_peaks,
                 aa_dict_nun)
        
    # print ("All info from peaks for {} rates are in file: {}".format(ccr_name, new_list), file=RaportBox)
    RaportBox.write("All info about zero CCR rates are in file: {}\n".format(new_list))
    return


def deepResidueAnalisys(writer,one_peak:CResidue,peak_number:int,experiment:CSpectrum,
                 error_residue_dict:dict,
                 error_experiment_dict:dict,
                 error_aa_dict:dict):
    """ For particular peak we write all information about it and append basics info inf dicts
    
    """
    s_dim = experiment.n_dim
    CCR_name =  experiment.CCR_name+Additional_text(experiment)
    noise_level = experiment.noise
    peak_descr = str(one_peak.aa_number)+one_peak.aa_name
    list_of_row = []
    one_row = {}
    error_experiment_dict[CCR_name]+=1
    if one_peak.aa_name in error_aa_dict:
        error_aa_dict[one_peak.aa_name]+=1
    else:
        error_aa_dict[one_peak.aa_name]=1
    one_row["CCR name"] = CCR_name
    one_row["AA"] = peak_descr
    if one_peak.is_peak == True and one_peak.is_ccr_rate == True:
        if experiment.CCR_pos == 0:
            for i in range(s_dim):
                one_row['w{}'.format(i+1)] = '{:.3f}'.format(one_peak.peak_pos[i])
            one_row["intensity in auto"] = one_peak.peak_intens[0]
            one_row["intensity in cross"] = one_peak.peak_intens[1]
            one_row["Noise in auto"] = noise_level[0]
            one_row["Noise in cross"] = noise_level[1]

        if one_peak.ccrrate_calculation_error == False:
            one_row["CCR rate"] = '{:.4f}'.format(one_peak.ccr_rate)
        else:
            one_row["Comments"] = "peak intens error"

    elif one_peak.is_peak == False and one_peak.is_ccr_rate == True: 
        if one_peak.aa_name == "G":
            one_row["Comments"] = "glycine - no CCR rate"
        else:
            one_row["CCR rate"] = '{:.4f}'.format(one_peak.ccr_rate)
            one_row["Comments"] = "Info from other peak"
    elif one_peak.is_overlap == True:
        one_row["Comments"] = "This peak overlap with: {}".format(one_peak.overlap_peaks)
    elif one_peak.is_peak == True and one_peak.is_ccr_rate == False:
        for i in range(s_dim):
            one_row['w{}'.format(i+1)] = '{:.3f}'.format(one_peak.peak_pos[i])
    else:
        one_row["Comments"] = "No visible peak"
    
    if one_peak.ccrrate_error_value != 0.0:
        one_row['Uncertainty'] = '{:.4f}'.format(one_peak.ccrrate_error_value)            # ERROR GAMMA INFO, 9.11.2023
    

    if one_peak.is_gamma_calc and one_peak.is_ccr_rate:
        one_row['Reference Gamma'] = "{:.4f}".format(one_peak.gamma_ref)
        one_row['Ix Theor'] = "{:.0f}".format(one_peak.Ix_theor)
        one_row['Ix/Ia'] = "{:.4f}".format(one_peak.Ix_theor_Ia_ratio)
        one_row['Ix/Ia without NS'] = "{:.4f}".format(one_peak.Ix_theor_Ia_ratio_without_NS)
    
    writer.writerow(one_row)

    if experiment.CCR_pos == -1:
        next_peak = experiment.peak[peak_number+1]
        one_row = {}
        peak_descr = str(next_peak.aa_number)+next_peak.aa_name
        one_row["CCR name"] = CCR_name
        one_row["AA"] = peak_descr
        for i in range(s_dim):
            one_row['w{}'.format(i+1)] = '{:.3f}'.format(next_peak.peak_pos[i])
        one_row["intensity in auto"] = next_peak.peak_intens[0]
        one_row["intensity in cross"] = next_peak.peak_intens[1]
        one_row["Noise in auto"] = noise_level[0]
        one_row["Noise in cross"] = noise_level[1]
        writer.writerow(one_row)

        if peak_descr in error_residue_dict:
            error_residue_dict[peak_descr].append(deepcopy(CCR_name+" (-1 position)"))
        else:
            error_residue_dict[peak_descr] = [CCR_name+" (-1 position)"]
    elif experiment.CCR_pos == 0:
        if peak_descr in error_residue_dict:
            error_residue_dict[peak_descr].append(deepcopy(CCR_name))
        else:
            error_residue_dict[peak_descr] = [CCR_name]


def deepVisibilityAnalisis(one_peak:CResidue,experiment:CSpectrum,invisible_cross_peaks):
    """Checking if peak if visible in cross version of experiment
    """
    CCR_name =  experiment.CCR_name+Additional_text(experiment)
    noise_level = experiment.noise
    peak_descr = str(one_peak.aa_number)+one_peak.aa_name
    if one_peak.is_peak:
        if abs(one_peak.peak_intens[1])<noise_level[1]*10:
            invisible_cross_peaks[CCR_name][0]+=1
            RaportBox.write("{} - {}, noise: {} intens: {} \n".format(invisible_cross_peaks[CCR_name][0],peak_descr,noise_level[1],one_peak.peak_intens[1]))
        if abs(one_peak.peak_intens[1])<noise_level[1]:
            invisible_cross_peaks[CCR_name][1]+=1

def deepExpAnalisys(writer_to_file_by_dict:csv,
                 error_residue_dict:dict,
                 error_experiment_dict:dict,
                 error_aa_dict:dict,
                 invisible_cross_peaks:dict,
                 aa_dict_nun:dict):
    """ 
    1) For peak in error_residue_dict are write down inforation about what residues in which CCR experiment they are
    2) Per CCR experiment are write down information about how many peaks are error and if its intensity if smaller than 10 time noise val or noise val
    3)
    """
    for residue in sorted(error_residue_dict.keys()):
        output_line = [residue]
        for ccr_exp in error_residue_dict[residue]:
            output_line.append(deepcopy(ccr_exp))
        writer_to_file_by_dict.writerow(output_line)
    writer_to_file_by_dict.writerow("")
    writer_to_file_by_dict.writerow(["Number of peaks in cross where:","CCR zero error","intesity < 10*noise level","intesity < noise level"])
    for experiment in Experiments:
        CCR_name =  experiment.CCR_name+Additional_text(experiment)
        writer_to_file_by_dict.writerow([CCR_name,error_experiment_dict[CCR_name],invisible_cross_peaks[CCR_name][0],invisible_cross_peaks[CCR_name][1]])
    writer_to_file_by_dict.writerow("")
    writer_to_file_by_dict.writerow(["AA name","CCR zero error"])
    for aa in aminoacids:
        if Res3to1(aa) in error_aa_dict:
            # print_raport(f"{aa},{error_aa_dict[Res3to1(aa)]},{aa_dict_nun[Res3to1(aa)]}\n")
            percent = error_aa_dict[Res3to1(aa)]/aa_dict_nun[Res3to1(aa)]*100
            writer_to_file_by_dict.writerow([aa,error_aa_dict[Res3to1(aa)],f"{percent:.1f} %"])



def Write_ErrorCCRrates_Analisis_CSV():
    new_list = "ErrorCCRrates_Analisis.csv".format(RaportDir)
    
    s_dim = Experiments[0].n_dim
    error_residue_dict={}
    error_experiment_dict={}
    error_aa_dict = {}
    invisible_cross_peaks = {}
    aa_dict_nun = {}

    with open(new_list, mode='w', newline='') as csv_file:
        headers = ['CCR name','AA']
        for i in range(1,s_dim+1):
            headers.append(deepcopy('w{}'.format(i)))
        headers.extend(deepcopy(['intensity in auto','Noise in auto','intensity in cross','Noise in cross','CCR rate','Uncertainty','Comments','Reference Gamma', 'Ix Theor', 'Ix/Ia', 'Ix/Ia without NS']))
        writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
        writer.writeheader()
        for experiment in Experiments:
            CCR_name =  experiment.CCR_name+Additional_text(experiment)
            if experiment.is_peaklist == True:
                writer.writerow({})
                invisible_cross_peaks[CCR_name]=[0,0]
                error_experiment_dict[CCR_name]=0
                for peak_number, one_peak in enumerate(experiment.peak):
                    if one_peak.aa_name in aa_dict_nun:
                        aa_dict_nun[one_peak.aa_name]+=1
                    else:
                        aa_dict_nun[one_peak.aa_name]=1
                    deepVisibilityAnalisis(one_peak,experiment,invisible_cross_peaks)

                    if one_peak.is_ccr_rate and one_peak.aa_name != "G" and one_peak.gamma_ref > one_peak.ccr_rate*5 and one_peak.is_gamma_calc:
                            deepResidueAnalisys(writer,one_peak,peak_number,experiment,
                                                    error_residue_dict,
                                                    error_experiment_dict,
                                                    error_aa_dict)

        writer_row = csv.writer(csv_file,delimiter=",")
        writer_row.writerow("")
        
        deepExpAnalisys(writer_row,
                 error_residue_dict,
                 error_experiment_dict,
                 error_aa_dict,
                 invisible_cross_peaks,
                 aa_dict_nun)
        
        
    # print ("All info from peaks for {} rates are in file: {}".format(ccr_name, new_list), file=RaportBox)
    RaportBox.write("All info about zero CCR rates are in file: {}\n".format(new_list))
    return

def Write_diffCCRrates_Analisis_CSV():
    new_list = "diffCCRrates_Analisis.csv".format(RaportDir)
    
    s_dim = Experiments[0].n_dim
    error_residue_dict={}
    error_experiment_dict={}
    error_aa_dict = {}
    invisible_cross_peaks = {}
    aa_dict_nun = {}

    with open(new_list, mode='w', newline='') as csv_file:
        headers = ['CCR name','AA']
        for i in range(1,s_dim+1):
            headers.append(deepcopy('w{}'.format(i)))
        headers.extend(deepcopy(['intensity in auto','Noise in auto','intensity in cross','Noise in cross','CCR rate','Uncertainty','Comments','Reference Gamma', 'Ix Theor', 'Ix/Ia', 'Ix/Ia without NS']))
        writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
        writer.writeheader()
        for experiment in Experiments:
            CCR_name =  experiment.CCR_name+Additional_text(experiment)
            if experiment.is_peaklist == True:
                writer.writerow({})
                invisible_cross_peaks[CCR_name]=[0,0]
                error_experiment_dict[CCR_name]=0
                for peak_number, one_peak in enumerate(experiment.peak):
                    if one_peak.aa_name in aa_dict_nun:
                        aa_dict_nun[one_peak.aa_name]+=1
                    else:
                        aa_dict_nun[one_peak.aa_name]=1
                    deepVisibilityAnalisis(one_peak,experiment,invisible_cross_peaks)

                    if one_peak.is_ccr_rate and one_peak.aa_name != "G" and one_peak.is_gamma_calc:
                        if abs(one_peak.gamma_ref - one_peak.ccr_rate)>=3:
                            deepResidueAnalisys(writer,one_peak,peak_number,experiment,
                                                    error_residue_dict,
                                                    error_experiment_dict,
                                                    error_aa_dict)

        writer_row = csv.writer(csv_file,delimiter=",")
        writer_row.writerow("")
        
        deepExpAnalisys(writer_row,
                 error_residue_dict,
                 error_experiment_dict,
                 error_aa_dict,
                 invisible_cross_peaks,
                 aa_dict_nun)
        
        
    # print ("All info from peaks for {} rates are in file: {}".format(ccr_name, new_list), file=RaportBox)
    RaportBox.write("All info about zero CCR rates are in file: {}\n".format(new_list))
    return


"""                      MAIN PROGRAM                      """

if __name__ == "__main__":
    RaportBox = open(RaportBoxFile,"w")
    print_raport("\n=== Reading input files ===")
    Experiments, ToCompereDict = ReadExpSet(exp_file_name)
    Sequence_List, NumberOfAAPerResidue = ReadSequnce(seq_file_name)
    if refgammaFlag:
        gamma_ref_dict = Read_Reference_Gamma(gamma_cal_file_name)
    else: gamma_ref_dict = {}
    
    for experiment in Experiments:
        Residue_List = deepcopy(Sequence_List)
        Residue_List, experiment.is_peaklist, experiment.noise[0], experiment.noise[1] = Read_peaklist(file_director, experiment.auto_name, experiment.cross_name, experiment.n_dim, Residue_List, experiment.is_peaklist, points_mode=True)
        Residue_List, experiment.is_peaklist = Read_peaklist(file_director, experiment.auto_name, experiment.cross_name, experiment.n_dim, Residue_List, experiment.is_peaklist, points_mode=False)
        experiment.peak = Residue_List

        PeakUncertaintyFlag1 = Read_peak_uncertainty(file_director, experiment.auto_name, experiment.n_dim, experiment.peak, version=0)
        PeakUncertaintyFlag2 = Read_peak_uncertainty(file_director, experiment.cross_name, experiment.n_dim, experiment.peak, version=1)
        if PeakUncertaintyFlag1 and PeakUncertaintyFlag2:
            experiment.peak_uncertainty = True

    print_raport ("=== Reading input files finished ===\n\n=== Calculating CCR rates starting ===")
    for experiment in Experiments:
        print_raport("\nExperiment number: {}\ndim number:{}\nCCr posistion:{}\nnumber of scanes:{}\ntc time: {}\nnoise: {}\nother: {}".format(experiment.CCR_name, experiment.n_dim, 
                                                                                                                        experiment.CCR_pos, experiment.ns, 
                                                                                                                        experiment.tc_vol, experiment.noise,
                                                                                                                        experiment.other))
        if experiment.is_peaklist == True:
            Add_text = Additional_text(experiment)
            CheckOverlap(experiment.peak)
            CalcCCRRate(experiment.peak,experiment.CCR_pos,experiment.ns,experiment.tc_vol,experiment.CCR_name,experiment.noise,experiment.peak_uncertainty)
            print_raport("Calculation of CCR rates for {}{} is ready".format(experiment.CCR_name,Add_text))
            if refgammaFlag:
                if experiment.CCR_name in gamma_ref_dict:
                    Add_ref_gamma(gamma_ref_dict["seq_num"],gamma_ref_dict[experiment.CCR_name],experiment.peak)
                    calc_theor_Ix(experiment.peak, experiment.CCR_pos,experiment.tc_vol,experiment.ns,experiment.CCR_name)
            WriteCCRRate_small(experiment.peak, experiment.CCR_name,add=Add_text)
            WriteCCRRate_all_info(experiment.peak, experiment.CCR_name, experiment.n_dim, experiment.auto_name, experiment.cross_name, experiment.tc_vol,experiment.ns,add=Add_text)
            WriteCCRRate_all_info_CSV(experiment.peak, experiment.CCR_name, experiment.n_dim, experiment.auto_name, experiment.cross_name, experiment.tc_vol,add=Add_text)
    print_raport ("\n=== Calculating CCR rates finished ===")
    Write_ALL_CCRRate_CSV()
    if refgammaFlag:
        FlagNoComperison = False
        print_raport ("\n=== Compering with reference values of CCR rates ===\n")
        for CCR_type in ToCompereDict:
            if len(ToCompereDict[CCR_type])>1:
                PeaksTable = []
                DescripTable = []
                for exp_number in ToCompereDict[CCR_type]:
                    PeaksTable.append(deepcopy(Experiments[exp_number].peak))
                    DescripTable.append(deepcopy(Additional_text(Experiments[exp_number])))
                    # plot_gamma_gamma(Experiments[exp_number].peak,CCR_type,
                    #                  add=Additional_text(Experiments[exp_number]),
                    #                  add2=Additional_text(Experiments[exp_number]))
                    # plot_error_histogram(Experiments[exp_number].peak,CCR_type,
                    #                  add=Additional_text(Experiments[exp_number]),
                    #                  add2=Additional_text(Experiments[exp_number]))
                # plot_gamma_diff(PeaksTable, CCR_type, DescripTable)
                plot_gamma_gamma_and_diff_together(PeaksTable, CCR_type, DescripTable)
                plot_cross_intens_theor_exp_together(PeaksTable, CCR_type, DescripTable)
                plot_error_histogram_together(PeaksTable, CCR_type, DescripTable)
                FlagNoComperison = False
            else:
                plot_gamma_gamma(Experiments[ToCompereDict[CCR_type][0]].peak,CCR_type,add=Additional_text(Experiments[ToCompereDict[CCR_type][0]]))
                plot_error_histogram(Experiments[ToCompereDict[CCR_type][0]].peak,CCR_type,add=Additional_text(Experiments[ToCompereDict[CCR_type][0]]))
                FlagNoComperison = True
        if FlagNoComperison:
            PeaksTable = []
            DescripTable = []
            CCRTable = []
            for exp in Experiments:
                PeaksTable.append(deepcopy(exp.peak))
                DescripTable.append(deepcopy(Additional_text(exp)))
                CCRTable.append(deepcopy(exp.CCR_name))
            plot_gamma_gamma_all_together(PeaksTable, CCRTable, DescripTable,publication_style=PublicationFlag)
        Write_zeroCCRrates_Analisis_CSV()
        Write_ErrorCCRrates_Analisis_CSV()
        Write_diffCCRrates_Analisis_CSV()
        print_raport ("=== Compering with reference values of CCR rates finished ===\n")
    else:
        print_raport ("\n=== Compering experiments with diffrent number of NUS points ===\n")
        for CCR_type in ToCompereDict:
            if len(ToCompereDict[CCR_type])>1:
                max_number = 0.0
                min_number = 10000000.0
                min_exp = -1
                max_exp = -1 
                for exp_number in ToCompereDict[CCR_type]:
                    nus_number = Additional_text(Experiments[exp_number])[1:-3]
                    if nus_number.isdigit() == True:
                        nus_number = float(nus_number)
                        # if nus_number.is_integer:
                        if nus_number > max_number:
                            max_number = nus_number
                            max_exp = exp_number
                        if nus_number < min_number:
                            min_number = nus_number
                            min_exp = exp_number
                if max_number.is_integer:
                    plot_gamma_exp_vs_epx(Experiments[min_exp].peak,Experiments[max_exp].peak, CCR_type,str(int(min_number)),str(int(max_number)))
                else:
                    plot_gamma_exp_vs_epx(Experiments[min_exp].peak,Experiments[max_exp].peak, CCR_type,str(min_number),str(max_number))
        print_raport ("\n=== Compering experiments with diffrent number of NUS points finished ===\n")
    RaportBox.close()





    # TO DO
    # dodać predefiniowane eksperymenty
    # dodac ograniczenie dla numerów reszt dla których liczymy 
    # dodac dodawanie daty i adnotacji do wykresow 