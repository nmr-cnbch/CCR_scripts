#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 27 11:18 2022

@author: paulina
"""
""" Program do czytania list pików i obliczania stałych CCR"""
    
    
from math import atanh
from copy import deepcopy
import sys

if len(sys.argv)==1:
    sys.exit("""
    For run script type in command line:
        python3 calc_CCR_rate [file director]

    additionaly you can add:
        --  """)

file_director = sys.argv[1]
if file_director[-1] != "/":
    print (file_director)
    file_director += "/"


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
        self.tc_vol = -1        # time of ccr evolution
        self.Hroi = [0,0]       # downfield and upfield of direct dimention

        self.peak = []            #informations about peaks from CResidue class
            
                    
class CResidue:
    def __init__(self):
        self.aa_name = "None"
        self.aa_number = 0              # the number of the amino acid to which this peak corresponds
        self.is_peak = False       # information about the presence of a peak


        self.peak_pos = []            # chemical shifts for all nuclei of peak, length depends on dimentionality
        self.peak_pos_points = []
        self.peak_intens = [0,0]      # peak hight in auto [0] and cross [1] version
        self.descript = ""              # description of peak which is in first column in Sparky-like peak list
        # self.peak_seq_pos_a = []

        self.is_overlap = False       # information about the presents of overlaping of a peak: maybe, yes, no   
        self.overlap_peaks = []       #      
        self.is_ccr_rate = False
        self.ccr_rate = -99999


        
        




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
        sum_of_squares += (peak1[i]-peak2[i])**2
    dis =  sum_of_squares**(1/len(peak1))
    return dis

def ReadExpSet(file_director):               # wczytywanie danych z pliku experiments_set.txt do klasy CSpectrum
    experiments = []
    
    with open("{}experiments_set.txt".format(file_director), "r") as exp_set:
        print ("otwarte")
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
        print ("expset_lines", expset_lines)
        if commentFlag == False:
            expset_lines.append(deepcopy(len(lines)))

        for i in range(0,len(expset_lines)-1):
            # print ("Lines between ", expset_lines[i], expset_lines[i+1])
            one_experinet = CSpectrum()
            one_experinet.name=lines[expset_lines[i]][8:-1]
            for line in range(expset_lines[i],expset_lines[i+1]):
                if lines[line] != "\n":
                    items=lines[line].split()
                    if "type_of_CCR" in lines[line]:
                        one_experinet.CCR_name=items[1][1:-1]
                        print ("CCR_NAME", items[1][1:-1])
                    if "dir_auto" in lines[line]:
                        auto_name_dir = items[1]+"_"+one_experinet.CCR_name+"_a"
                        print ("auto_name_dir", auto_name_dir)
                        one_experinet.auto_name = auto_name_dir
                    if "dir_cross" in lines[line]:
                        cross_name_dir = items[1]+"_"+one_experinet.CCR_name+"_x"
                        print ("cross_name_dir", cross_name_dir)
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
            experiments.append(deepcopy(one_experinet))
            print ("numer_eksperymentu", i, "objekty_klasy", one_experinet.CCR_name, one_experinet.n_dim, one_experinet.nucl_name, one_experinet.nucl_pos, one_experinet.n_angle, one_experinet.angle_pos, one_experinet.angle, one_experinet.ns, one_experinet.tc_vol, one_experinet.Hroi)

    return experiments


def ReadSequnce(file_director):

    sequence = []
    with open("{}seq".format(file_director), "r") as input_file:
        linia_seq=input_file.readlines()
        FASTAFlag=False
        if ">" in linia_seq[0]:
            FASTAFlag=True
            print ("FASTAFlag")
            # FASTAseq=[]
            aa_no = 0
            for indexl,l_seq in enumerate(linia_seq[1:]):
                oneletter=list(str(l_seq))
                for aa in oneletter: 
                    if aa.isalpha():
                        seq = CResidue()
                        seq.aa_name = aa         
                        aa_no += 1
                        seq.aa_number = aa_no
                        sequence.append(deepcopy(seq))
                        # FASTAseq.append(deepcopy(seq))
            
            ## Printing sequence     
            # print ("seqence lenth: ", len(sequence))
            # print ("seqence: ")
            # for indexf, fastaseq in enumerate(sequence):
            #     # print (fastaseq.aa_name)            
            #     aa_rest=Res1to3(fastaseq.aa_name) # check whether seq file contains correct aa names 
            #     print ("%d%s"%(indexf+1, aa_rest))
    return sequence


def Read_peaklist(file_director, peak_list_name_auto, peak_list_name_cross, rel_pos, s_dim, p_list, points_mode=False):
    
    NameFlag = [False, False]
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
                f.close
            except FileNotFoundError:
                print ("There is no such file or directory:", peaklistfile)
    

    if NameFlag[0] == True and NameFlag[1] == True:
        with open(peak_list_names[0], 'r') as pl_a, open(peak_list_names[1], 'r') as pl_x:  
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
                    try:
                        aminoacids_number = int(aminoacids[s_dim-1][1:-2])
                    except:
                        aminoacids_number = int(aminoacids[0][1:-1])
                    for indexr, residue in enumerate(p_list):
                        if aminoacids_number == residue.aa_number:
                            residue.is_peak = True
                            residue.descript = description
                            if points_mode == True:  # for points mode
                                for i in range(1,s_dim+1):
                                    residue.peak_pos_points.append(deepcopy(int(item_a[i])))
                            else: 
                                for i in range(1,s_dim+1):
                                    residue.peak_pos.append(deepcopy(float(item_a[i])))
                    # Reading peak intensity
                            if len(item_a)>s_dim+1:
                                residue.peak_intens[0] = float(item_a[s_dim+1])
                                residue.peak_intens[1] = float(item_x[s_dim+1])
    else:
        print ("There is no such file or directory:", peak_list_names)
    return p_list



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
                            peak1.overlap_peaks.append(peak2.aa_number)
                            peak2.overlap_peaks.append(peak1.aa_number)
    return



def CalcCCRRate(peak_list,angle_position,number_scan,Tc):
    for indexp, peak in enumerate(peak_list):
        if not peak.is_overlap and peak.is_peak:
            peak_list[indexp+angle_position].is_ccr_rate = True
            ccr_rate_vol = -atanh((peak.peak_intens[1]*number_scan[0])/(peak.peak_intens[0]*number_scan[1]))/Tc
            peak_list[indexp+angle_position].ccr_rate = ccr_rate_vol
    return



def WriteCCRRate(peak_list, file_director, ccr_name, s_dim):
    new_list = "{0}{1}_CCRrate.list".format(file_director,ccr_name)
    print ("new_peak_list", new_list)
    max_lenth_discrip = 0
    max_lenth_intens = 0
    for p in peak_list:
        pd = str(p.aa_number)+p.aa_name
        if len(pd)>max_lenth_discrip:
            max_lenth_discrip=len(pd)
        if len(str(p.peak_intens[0]))>max_lenth_intens:
            max_lenth_intens=len(str(p.peak_intens[0]))
    with open(new_list, 'w') as listfile:
        print ("\tInformation of files", file=listfile) 
        for one_peak in peak_list:
            if one_peak.is_ccr_rate == True:
                peak_descr = str(one_peak.aa_number)+one_peak.aa_name
                print ("{:{sentence_len}}".format(peak_descr, sentence_len=max_lenth_discrip), end="\t", file=listfile)
                if one_peak.is_peak == True:
                    for i in range(s_dim):
                        print ("{:.3f}".format(one_peak.peak_pos[i]), end="\t", file=listfile)
                    print ("{:{sentence_len}}\t{:{sentence_len}}\t{:.4f}".format(one_peak.peak_intens[0], one_peak.peak_intens[1], one_peak.ccr_rate, sentence_len=max_lenth_intens), file=listfile)
                else: 
                    print ("{:{sentence_len}}\t{:{sentence_len}}\t{:.4f}".format("Info from other peak","",one_peak.ccr_rate,sentence_len=7*s_dim), file=listfile)
                
    return



print ("start")

Experiments = ReadExpSet(file_director)
Residue_List = ReadSequnce(file_director)
for experiment in Experiments:
    Read_peaklist(file_director, experiment.auto_name, experiment.cross_name, experiment.nucl_pos, experiment.n_dim, Residue_List, True)
    Read_peaklist(file_director, experiment.auto_name, experiment.cross_name, experiment.nucl_pos, experiment.n_dim, Residue_List, False)
    experiment.peak = Residue_List
    CheckOverlap(experiment.peak)
    if len(experiment.angle_pos)==1 or experiment.angle_pos[0] == experiment.angle_pos[1]:
        CalcCCRRate(experiment.peak,experiment.angle_pos[0],experiment.ns,experiment.tc_vol)
    WriteCCRRate(experiment.peak, file_director, experiment.CCR_name, experiment.n_dim)


