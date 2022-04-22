#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 27 11:18 2022

@author: paulina
"""
""" Program do czytania list pików i obliczania stałych CCR"""
    
    
from copy import deepcopy
import sys

if len(sys.argv)==1:
    sys.exit("""
    For run script type in command line:
        python3 calc_CCR_rate [peak list1 path] [peak list2 path] [peak list3 path]...[peak listN path] --dim []
                        
        --dim [] - dimentionality of peak list

    additionaly you can add:
        --comp2list [num] [num]\t- to additionaly compare two of peak list from above
        --out_dir [path]\t- path to direction where output file will be put
        --fl [path]\t- path to file with list of relativ path to peaklists - every set of file must separate by ------\n""")


file_director = sys.argv[1]



class CSpectrum:
    def __init__(self):
        self.CCR_name = ''       # name of CCR: CCR_1, CCR_2
        self.auto_name = ''       # name of file with auto version
        self.cross_name = ''       # name of file with cross version
        self.n_dim = 0            # number of dimentions
        self.nucl_name = []            # nuclei of all peaks: H, N, C, CA, CB, HA, HB
        self.nucl_pos = []             # nuclei position of all peaks: -2, -1, 0, 1, 2 
        
        self.n_angle = 0          # number of measured angle  
        self.angle = []           # angle: phi, psi
        self.angle_pos = []       # angle position of all peaks: -2, -1, 0, 1, 2
        self.ns = [0,0]          # number of scans in auto [0] and cross [1] version
        self.tc_vol = -1        # time of ccr evolution
        self.Hroi = [0,0]       # downfield and upfield of direct dimention

        self.peak = []            #informations about peaks from CPeak class
            
                    
class CPeak:
    def __init__(self):
        self.is_peak = False       # information about the presence of a peak
        self.peak_pos = []            # chemical shifts for all nuclei of peak, length depends on dimentionality
        self.peak_intens = [0,0]      # peak hight in auto [0] and cross [1] version
        self.descript = ""              # description of peak which is in first column in Sparky-like peak list
        self.aa_number = 0              # the number of the amino acid to which this peak corresponds

        self.is_overlap = ''       # information about the presents of overlaping of a peak: maybe, yes, no   
        self.overlap_peaks = []       #      




class CSequence:
    def __init__(self):
        self.aa_name="None"
        self.aa_num=-1
        self.ccr_rate=-99999




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









def ReadExpSet(file_director):               # wczytywanie danych z pliku experiments_set.txt do klasy CSpectrum
    experiments = []
    
    with open("{}/experiments_set.txt".format(file_director), "r") as exp_set:
        print ("otwarte")
        lines = exp_set.readlines()
        expset_lines=[]
        for indexl, line in enumerate(lines):
            if "type_of_CCR" in line:
                expset_lines.append(deepcopy(indexl))
        print ("expset_lines", expset_lines)
        expset_lines.append(deepcopy(len(lines)))

        for i in range(0,len(expset_lines)-1):
            # print ("Lines between ", expset_lines[i], expset_lines[i+1])
            one_experinet=CSpectrum()
            one_experinet.name=lines[expset_lines[i]][8:-1]
            for line in range(expset_lines[i],expset_lines[i+1]):
                if len(lines[line])>0:
                    items=lines[line].split()
                    if "type_of_CCR" in lines[line]:
                            one_experinet.CCR_name=items[1]
                    # if "dir_auto" in lines[line]:
                    # 
                    #         one_experinet.auto_name=items[1]
                    # if "dir_cross" in lines[line]:
                    # 
                    #         one_experinet.cross_name=items[1]
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
            FASTAseq=[]
            aa_no = 0
            for indexl,l_seq in enumerate(linia_seq):
                if 0<indexl<len(linia_seq):
                    oneletter=list(str(l_seq))
                    for aa in oneletter: 
                        if aa.isalpha():
                            seqq = CSequence
                            seqq.aa_name = aa
                            print ("pierwsza",seqq.aa_name)            
                            aa_no += 1
                            # seq.aa_num = aa_no
                            sequence.append(seqq)
                            print ("druga",sequence[aa_no-1].aa_name, len(sequence))
                            FASTAseq.append(deepcopy(seqq))
                    # for i in range(len(oneletter)):
                    #     print ("druga",sequence[i].aa_name, len(sequence))
                    #     print ("trzecia",FASTAseq[i].aa_name, len(FASTAseq))
                            
            print ("seqence lenth: ", len(sequence))
            print ("seqence: ")
            for indexf, fastaseq in enumerate(sequence):
                print (fastaseq.aa_name)            
                aa_rest=Res1to3(fastaseq.aa_name) # check whether seq file contains correct aa names 
                print ("%d%s"%(indexf+1, aa_rest))
    return sequence, aa_no




def Read_peaklist(file_director, peak_list, s_dim):
    peaklistdir = peak_list.auto_name
    with open("{}{}".format(file_director,peaklistdir), 'r') as pl:  
        # print (peak_list)
        p_lines = pl.readlines()
        p_list = []
        aa_max_number = 1
        for indexl, line in enumerate(p_lines):
            if indexl > 1 :
                p_pos = CPeak()
                item = line.split()
                p_pos.descript = item[0]
                aa = p_pos.descript.split("-")
                try:
                    p_pos.aa_number = int(aa[s_dim-1][1:-2])
                except:
                    p_pos.aa_number = int(aa[0][1:-1])
                if p_pos.aa_number>aa_max_number:
                    aa_max_number=p_pos.aa_number
                if item[s_dim].isdigit():
                    # p_pos.peak_pos.append(deepcopy(int(item[s_dim])))
                    for i in range(1,s_dim+1):
                        p_pos.peak_pos.append(deepcopy(int(item[i])))
                else: 
                    # p_pos.peak_pos.append(deepcopy(float(item[s_dim])))
                    for i in range(1,s_dim+1):
                        p_pos.peak_pos.append(deepcopy(float(item[i])))
                if len(item)>s_dim+1:
                    p_pos.peak_intens = float(item[s_dim+1])
                p_list.append(deepcopy(p_pos))
    return p_list, aa_max_number




def ReadInputFile():
    return


def CheckOverlap():
    return


def CalcCCRRate():
    return


def WriteCCRRate():
    return


print ("start")
Experiments = ReadExpSet(file_director)
AASequence = ReadSequnce(file_director)