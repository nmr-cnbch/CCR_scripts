#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 27 11:18 2022

@author: paulina
"""
""" Program do czytania list pików i obliczania stałych CCR"""
    
    
from math import atanh
from copy import deepcopy
import math
import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy import stats
import os

if len(sys.argv)==1:
    sys.exit("""
    The script to calculate CCR constants. In director has to contains:
    - peak lists (with peak position in points and peak height)
    - experiment_set
    - sequence in FASTA format

    For run script type in command line:
        python3 calc_CCR_rate.py [file director]

    additionaly you can add:
        --seq  - if name of file with aminoacid sequence is not 'seq', add this with name of file
        --refgamma  - if you have file with reference values of CCR rates, add this with name of file
                      (file must be .csv, columns name should be: AA, psi_angle (or/and phi_angle), CCR name
        --expset    - if you want use experiments setup file with diffrent filename than "experiments_set.txt" (structure of file must be the same as orginal file) """)


"""Commantline reading"""

file_director = sys.argv[1]
if file_director[-1] != "/":
    file_director += "/"
    print ("\nFile director",file_director)

if "--seq" in sys.argv:    
    i = sys.argv.index("--seq")
    seq_file_name = sys.argv[i+1] 
else: seq_file_name = "{}seq".format(file_director)
print("File with aminoacid sequence: {}".format(seq_file_name)) 

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
# 




"""Classes"""

class CSpectrum:
    def __init__(self):
        self.CCR_name = ""       # name of CCR: CCR_1, CCR_2
        self.auto_name = ""       # name of file with auto version
        self.cross_name = ""       # name of file with cross version
        self.n_dim = 0            # number of dimentions
        self.nucl_name = []            # nuclei of all peaks: H, N, C, CA, CB, HA, HB
        self.nucl_pos = []             # nuclei position of all peaks: -2, -1, 0, 1, 2 
        
        self.n_angle = 0          # number of measured angle  
        self.angle = []           # angle: phi, psi
        self.angle_pos = []       # angle position of all peaks: -2, -1, 0, 1, 2
        self.ns = [0,0]          # number of scans in auto [0] and cross [1] version
        self.tc_vol = -1.0        # time of ccr evolution
        self.Hroi = [0.0,0.0]       # downfield and upfield of direct dimention

        self.peak = []            #informations about peaks from CResidue class
        self.is_peaklist = False
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
        self.is_overlap = False       # information about the presents of overlaping of a peak: maybe, yes, no   
        self.overlap_peaks = {}       #     

        self.is_ccr_rate = False
        self.ccr_rate = -99999.0
        self.ccrrate_error = False

        self.is_gamma_calc = False      # information about the presence of a reference gamma
        self.gamma_ref = -99999.0       # value of gamma (CCR rate) which will use as reference data for experimental data, value additional file with calculated previosly gammas
        


        
        

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
    RaportBox.write(text_to_write+"\n")



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
            experiments.append(deepcopy(one_experinet))
            print ("{} - {}\n\t reference exp: {}\n\t transfer exp: {}".format(one_experinet.CCR_name,one_experinet.other,one_experinet.auto_name,one_experinet.cross_name))
            # print ("Experiment number:{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}".format(i, one_experinet.CCR_name, one_experinet.n_dim, one_experinet.nucl_name, one_experinet.nucl_pos, one_experinet.n_angle, one_experinet.angle_pos, one_experinet.angle, one_experinet.ns, one_experinet.tc_vol, one_experinet.Hroi), file=RaportBox)
            RaportBox.write("\nExperiment number:{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n".format(i, one_experinet.CCR_name, one_experinet.n_dim, one_experinet.nucl_name, one_experinet.nucl_pos, one_experinet.n_angle, one_experinet.angle_pos, one_experinet.angle, one_experinet.ns, one_experinet.tc_vol, one_experinet.Hroi))
    return experiments, to_compere_dict


def ReadSequnce(seq_file_name):

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
                        # FASTAseq.append(deepcopy(seq))
            print_raport ("FASTA sequence: {}".format(FASTAseq))
    return sequence


def Read_peaklist(file_director, peak_list_name_auto, peak_list_name_cross, s_dim, p_list, is_there_peaklist, points_mode=False):
    
    NameFlag = [False, False]
    # NotFoundName = []
    if points_mode == False:
        listofnames = [".list","_new_ppm.list"]
    else:
        listofnames = ["_points.list","_new_points.list"]
    peak_list_basic_names = [peak_list_name_auto, peak_list_name_cross]
    peak_list_names = ["None","None"]
    for indexlv, list_verson in enumerate(peak_list_basic_names):
        for i, l in enumerate(listofnames):
            peaklistfile = str(file_director+list_verson+l)
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
        # print ("\n\nLISTA:", peak_list_names)
        with open(peak_list_names[0], 'r') as pl_a, open(peak_list_names[1], 'r') as pl_x:  
            is_there_peaklist = True
            p_lines_a = pl_a.readlines()
            p_lines_x = pl_x.readlines()
            for indexl, line in enumerate(p_lines_a):
                if indexl > 1 :
                    item_a = line.split()
                    item_x = p_lines_x[indexl].split()
                    # Reading description
                    if item_a[0] != item_x[0]:
                        print ("Error: Auto ({}) and cross ({}) peak list are not compatible:\n{} ? {}".format(peak_list_names[0], peak_list_names[1], item_a[0], item_x[0]))
                        sys.exit()
                    description = item_a[0]
                    aminoacids = description.split("-")
                    # print (aminoacids)
                    try:
                        aminoacids_number = int(aminoacids[s_dim-1][1:-2])
                    except:
                        aminoacids_number = int(aminoacids[0][1:-2])
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
                                residue.peak_intens[0] = float(item_a[s_dim+1])
                                residue.peak_intens[1] = float(item_x[s_dim+1])
    else:
        print_raport ("Missing file for {} and {}! {}\n".format(peak_list_name_auto, peak_list_name_cross, peak_list_names))
        is_there_peaklist = False
    return p_list, is_there_peaklist

def Read_Reference_Gamma(gamma_cal_file_name):
    with open(gamma_cal_file_name, newline='') as gc_file:  
        reader = csv.DictReader(gc_file, delimiter="\t")
        headers = next(reader)
        gamma_file_dict = {}
        for item in headers:
            gamma_file_dict[item]=[headers[item]]
        for col in reader:
            for item in headers:
                # print (item)
                gamma_file_dict[item].append(col[item])
    
    try:
        Res3to1(gamma_file_dict["AA"][0][-3:])
        if gamma_file_dict["AA"][0][:-3].isnumeric():
            gamma_file_dict["seq_num"]=[gamma_file_dict["AA"][0][:-3]]
            for i in range(1,len(gamma_file_dict["AA"])):
                # print ("gamma_file_dict: AA - {}, NUM - {}".format(gamma_file_dict["AA"][i],gamma_file_dict["AA"][i][:-3]))
                gamma_file_dict["seq_num"].append(deepcopy(gamma_file_dict["AA"][i][:-3]))
    except:
        if gamma_file_dict["AA"][0][:-1].isnumeric():
            gamma_file_dict["seq_num"]=[gamma_file_dict["AA"][0][:-1]]
            for i in range(gamma_file_dict["AA"]):
                gamma_file_dict["seq_num"].append(deepcopy(gamma_file_dict["AA"][i][:-1]))

    return gamma_file_dict


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



def CalcCCRRate(peak_list,angle_position,number_scan,Tc,CCR_name):
    for indexp, peak in enumerate(peak_list):
        if not peak.is_overlap and peak.is_peak:
            peak_list[indexp+angle_position].is_ccr_rate = True
            try:
                # print (peak.peak_intens[1], number_scan[0], peak.peak_intens[0],number_scan[1], "---->", (peak.peak_intens[1]*number_scan[0])/(peak.peak_intens[0]*number_scan[1]))
                if CCR_name == "CCR_5" or CCR_name == "CCR_6":
                    ccr_rate_vol = -atanh((peak.peak_intens[1]*number_scan[0])/(peak.peak_intens[0]*number_scan[1]))/Tc
                else:
                    ccr_rate_vol = atanh((peak.peak_intens[1]*number_scan[0])/(peak.peak_intens[0]*number_scan[1]))/Tc
            except:
                ccr_rate_vol = (peak.peak_intens[1]*number_scan[0])/(peak.peak_intens[0]*number_scan[1]) #"atanh(x) - x shoudl be beetween -1 and 1"
                peak.ccrrate_error = True
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
            if int(one_peak.aa_number) == int(r_seq_num):
                one_peak.is_gamma_calc = True
                one_peak.gamma_ref = float(ref_CCR[indexr])
                # print ("gamma_ref",one_peak.gamma_ref)
    return



def LRegression_expresion(x,y):
    slope, intercept, r, p, std_err = stats.linregress(x, y)
    r2 = r**2
    lr_exp = '{:.3f}*x + {:.3f}'.format(slope,intercept)
    lr_exp_y_list = []
    for i in x:
        lr_exp_y_list.append(deepcopy(slope*i+intercept))
    return lr_exp_y_list, lr_exp, r2



def plot_gamma_gamma(peaks,ccr_name,transparent_plot=False,add=""):

    min_max_value = [+100.0,-100.0]         #minimal and maximal value of CCR rate or reference gamma - it is necessary to make plot  
    gamma_calculated = []
    gamma_experimental = []
    for one_peak in peaks:
        if one_peak.is_ccr_rate and one_peak.is_gamma_calc:
            gamma_calculated.append(deepcopy(one_peak.gamma_ref))
            gamma_experimental.append(deepcopy(one_peak.ccr_rate))
            if one_peak.gamma_ref<min_max_value[0]:
                min_max_value[0]=one_peak.gamma_ref
            if one_peak.ccr_rate<min_max_value[0]:
                min_max_value[0]=one_peak.ccr_rate
            if one_peak.gamma_ref>min_max_value[1]:
                min_max_value[1]=one_peak.gamma_ref
            if one_peak.ccr_rate>min_max_value[1]:
                min_max_value[1]=one_peak.ccr_rate

    # print ("gamma_calculated:\t", gamma_calculated)
    # print ("gamma_experimental:\t", gamma_experimental)
    # print ("min_max_value:\t", min_max_value)
        
    plt.rcParams['font.size'] = '14'

    # Plotting both the curves simultaneously
    plt.axline([0,0],slope=1, linestyle=(0, (5, 5)), linewidth=1.5, color='darkgray', label='x=y')
    plt.scatter(gamma_calculated, gamma_experimental, s=10, color='#0066ffff', label='structure-predicted vs experimental')
    
    # Calculating Linear Regression
    # lr_exp_y_list,lr_expresion, coefition, = LRegression_expresion(gamma_calculated,gamma_experimental)
    # plt.plot(gamma_calculated, lr_exp_y_list, color='forestgreen', label='{}, $R^2$ = {:.3f}'.format(str(lr_expresion),coefition))
    # plt.plot(gamma_calculated, lr_exp_y_list, color='orange', label='{}, $R^2$ = {:.3f}'.format(str(lr_expresion),coefition))
    
    
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
    # plt.legend(fontsize="small")
    
    # To load the display window
    # plt.savefig("{}/wykresy2/{}_exp_vs_calc.png".format(Dir, tfn[0].name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
    plt.savefig("{}{}_exp_vs_calc.png".format(file_director, ccr_name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
    # plt.show()
    plt.clf()
    return





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
                    plt.annotate(one_peak.aa_number,(one_peak.gamma_ref,one_peak.ccr_rate),textcoords="offset points",xytext=(0,-12),ha='center')
                if one_peak.ccr_rate-one_peak.gamma_ref>5 and indexp==0:
                    plt.annotate(one_peak.aa_number,(one_peak.gamma_ref,one_peak.ccr_rate),textcoords="offset points",xytext=(0,10),ha='center')
                if one_peak.gamma_ref<min_max_value[0]:
                    min_max_value[0]=one_peak.gamma_ref
                if one_peak.ccr_rate<min_max_value[0]:
                    min_max_value[0]=one_peak.ccr_rate
                if one_peak.gamma_ref>min_max_value[1]:
                    min_max_value[1]=one_peak.gamma_ref
                if one_peak.ccr_rate>min_max_value[1]:
                    min_max_value[1]=one_peak.ccr_rate
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
    return




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
            if peaks_min[indexp].ccr_rate<min_max_value[0]:
                min_max_value[0]=peaks_min[indexp].ccr_rate
            if one_peak.ccr_rate<min_max_value[0]:
                min_max_value[0]=one_peak.ccr_rate
            if peaks_min[indexp].ccr_rate>min_max_value[1]:
                min_max_value[1]=peaks_min[indexp].ccr_rate
            if one_peak.ccr_rate>min_max_value[1]:
                min_max_value[1]=one_peak.ccr_rate
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
    return





def WriteCCRRate_small(peak_list, ccr_name, add=""):
    new_list = "{}{}{}.csv".format(RaportDir,CCRname2Ratename(ccr_name),add)
    with open(new_list, mode='w', newline='') as csv_file:
        headers = ['AA','CCR rate']
        writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
        for one_peak in peak_list:
            one_row = {}
            one_row["AA"] = one_peak.aa_number
            if one_peak.is_ccr_rate == True and one_peak.ccrrate_error == False and one_peak.aa_name != "G":
                one_row["CCR rate"] = '{:.4f}'.format(one_peak.ccr_rate)
            else:
                one_row["CCR rate"] = "nan"
            writer.writerow(one_row)
    # print ("CCR rate for CCR efect between {} rates are in file: {}".format(CCRname2Ratename(ccr_name), new_list), file=RaportBox)
    RaportBox.write("CCR rate for CCR efect between {} rates are in file: {}\n".format(CCRname2Ratename(ccr_name), new_list))
    return




def WriteCCRRate_all_info(peak_list, ccr_name, s_dim, autoname, crossname, Tctime, add=""):
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
        print ("{}\t{:{sentence_len1}}\t{:{sentence_len2}}\t{:{sentence_len2}}\t{}\n".format("AA","peak position","intensity in auto","intensity in cross","CCR rate",sentence_len1=7*s_dim, sentence_len2=max_lenth_intens), file=listfile)
        for one_peak in peak_list:
            peak_descr = str(one_peak.aa_number)+one_peak.aa_name
            print ("{:{sentence_len}}".format(peak_descr, sentence_len=max_lenth_discrip,), end="\t", file=listfile)
            if one_peak.is_peak == True and one_peak.is_ccr_rate == True:
                for i in range(s_dim):
                    print ("{:.3f}".format(one_peak.peak_pos[i]), end="\t", file=listfile)
                if one_peak.aa_name == "G":
                    print ("{}".format("glycine - no CCR rate"), file=listfile)
                elif one_peak.ccrrate_error == False:
                    print ("{:{sentence_len}}\t{:{sentence_len}}\t{:.4f}".format(one_peak.peak_intens[0], one_peak.peak_intens[1], one_peak.ccr_rate, sentence_len=max_lenth_intens), file=listfile)
                else:
                    print ("{:{sentence_len}}\t{:{sentence_len}}\tpeak intens error: Ix/Ia = {:.4f}".format(one_peak.peak_intens[0], one_peak.peak_intens[1], one_peak.ccr_rate, sentence_len=max_lenth_intens), file=listfile)
                    
            elif one_peak.is_peak == False and one_peak.is_ccr_rate == True: 
                if one_peak.aa_name == "G":
                    print ("{}".format("\tglycine - no CCR rate"), file=listfile)
                else:
                    print ("{:{sentence_len}}\t\t\t\t\t\t{:.4f}".format("Info from other peak",one_peak.ccr_rate, sentence_len=7*s_dim), file=listfile)
            elif one_peak.is_overlap == True:
                print ("This peak overlap with: {}".format(one_peak.overlap_peaks), file=listfile)
            elif one_peak.is_peak == True and one_peak.is_ccr_rate == False:
                for i in range(s_dim):
                    print ("{:.3f}".format(one_peak.peak_pos[i]), end="\t", file=listfile)
                print ("", file=listfile)
            else:
                print ("{}".format("No visible peak"), file=listfile)
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
        headers.extend(deepcopy(['intensity in auto','intensity in cross','CCR rate','Comments']))
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
                    if one_peak.ccrrate_error == False:
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
                if one_peak.is_ccr_rate == True and one_peak.ccrrate_error == False and one_peak.aa_name != "G":
                    one_row['{}{}'.format(experiment.CCR_name,Add_textt)] = '{:.4f}'.format(one_peak.ccr_rate)
                else:
                    one_row['{}{}'.format(experiment.CCR_name,Add_textt)] = "nan"
            writer.writerow(one_row)
    print ("All CCR rates are in:", new_list)
    return



"""                      MAIN PROGRAM                      """

if __name__ == "__main__":
    RaportBox = open(RaportBoxFile,"w")
    print_raport("\n=== Reading input files ===")
    Experiments, ToCompereDict = ReadExpSet(exp_file_name)
    Sequence_List = ReadSequnce(seq_file_name)
    if refgammaFlag:
        gamma_ref_dict = Read_Reference_Gamma(gamma_cal_file_name)
    else: gamma_ref_dict = {}
    
    for experiment in Experiments:
        Residue_List = deepcopy(Sequence_List)
        Residue_List, experiment.is_peaklist = Read_peaklist(file_director, experiment.auto_name, experiment.cross_name, experiment.n_dim, Residue_List, experiment.is_peaklist, points_mode=True)
        Residue_List, experiment.is_peaklist = Read_peaklist(file_director, experiment.auto_name, experiment.cross_name, experiment.n_dim, Residue_List, experiment.is_peaklist, points_mode=False)
        experiment.peak = Residue_List

    print_raport ("=== Reading input files finished ===\n\n=== Calculating CCR rates starting ===")
    for experiment in Experiments:
        if experiment.is_peaklist == True:
            Add_text = Additional_text(experiment)
            CheckOverlap(experiment.peak)
            if len(experiment.angle_pos)==1 or experiment.angle_pos[0] == experiment.angle_pos[1]:
                CalcCCRRate(experiment.peak,experiment.angle_pos[0],experiment.ns,experiment.tc_vol,experiment.CCR_name)
                print_raport("Calculation of CCR rates for {}{} is ready".format(experiment.CCR_name,Add_text))
            if refgammaFlag:
                if experiment.CCR_name in gamma_ref_dict:
                    Add_ref_gamma(gamma_ref_dict["seq_num"],gamma_ref_dict[experiment.CCR_name],experiment.peak)
            WriteCCRRate_small(experiment.peak, experiment.CCR_name,add=Add_text)
            WriteCCRRate_all_info(experiment.peak, experiment.CCR_name, experiment.n_dim, experiment.auto_name, experiment.cross_name, experiment.tc_vol,add=Add_text)
            WriteCCRRate_all_info_CSV(experiment.peak, experiment.CCR_name, experiment.n_dim, experiment.auto_name, experiment.cross_name, experiment.tc_vol,add=Add_text)
    print_raport ("\n=== Calculating CCR rates finished ===")
    Write_ALL_CCRRate_CSV()
    print_raport ("\n=== Compering experiments with diffrent number of NUs points ===\n")
    for CCR_type in ToCompereDict:
        if len(ToCompereDict[CCR_type])>1:
            max_nus_number = 0
            min_nus_number = 10000000
            min_nus_exp = -1
            max_nus_exp = -1 
            for exp_number in ToCompereDict[CCR_type]:
                nus_number = Additional_text(Experiments[exp_number])[1:-3]
                if nus_number.isdigit:
                    if int(nus_number) > int(max_nus_number):
                        max_nus_number = int(nus_number)
                        max_nus_exp = exp_number
                    if int(nus_number) < int(min_nus_number):
                        min_nus_number = int(nus_number)
                        min_nus_exp = exp_number
            plot_gamma_exp_vs_epx(Experiments[min_nus_exp].peak,Experiments[max_nus_exp].peak, CCR_type,str(min_nus_number),str(max_nus_number))
    print_raport ("\n=== Compering experiments with diffrent number of NUs points finished ===\n")
    if refgammaFlag:
        print_raport ("\n=== Compering with reference values of CCR rates ===\n")
        for CCR_type in ToCompereDict:
            if len(ToCompereDict[CCR_type])>1:
                PeaksTable = []
                DescripTable = []
                for exp_number in ToCompereDict[CCR_type]:
                    PeaksTable.append(deepcopy(Experiments[exp_number].peak))
                    DescripTable.append(deepcopy(Additional_text(Experiments[exp_number])))
                plot_gamma_diff(PeaksTable, CCR_type, DescripTable)
            else:
                plot_gamma_gamma(Experiments[ToCompereDict[CCR_type][0]].peak,CCR_type,add=Additional_text(Experiments[ToCompereDict[CCR_type][0]]))
        print_raport ("=== Compering with reference values of CCR rates finished ===\n")
    RaportBox.close()