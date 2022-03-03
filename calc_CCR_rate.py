#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 27 11:18 2022

@author: paulina
"""
""" Program do czytania list pików i obliczania stałych CCR"""
    
    
from copy import deepcopy
from msilib import sequence
import sys



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

        self.is_overlap = ''       # information about the presents of overlaping of a peak: maybe, yes, no   
        self.overlap_peaks = []       #      




class CSequence:
    def __init__(self):
        self.aa_name=''
        self.aa_num=-1
        self.ccr_rate=-99999


def ReadExpSet():               # wczytywanie danych z pliku experiments_set.txt do klasy CSpectrum
    experiments = []
    
    with open("./experiments_set.txt", "r") as exp_set:
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


def ReadSequnce():

    sequence = []
    seq=CSequence

    with open("./seq", "r") as input_file:
        linia_seq=input_file.readlines()
        FASTAFlag=False
        if ">" in linia_seq[0]:
            FASTAFlag=True
            #print "FASTAFlag"
            FASTAseq=[]
            for indexl,l_seq in enumerate(linia_seq):
                if indexl+1<len(linia_seq):
                    oneletter=list(str(linia_seq[indexl+1]))
                    for indexo, aa in enumerate(oneletter): 
                        if oneletter[indexo].isalpha():
                            FASTAseq.append(deepcopy(oneletter[indexo]))
            print "seqence: ",
            for indexf, fastaseq in enumerate(FASTAseq):            
                iiii=Res1ToNr(fastaseq) # check whether seq file contains correct aa names 
                seq.res_type.append(deepcopy(Res1to3(fastaseq)))
                seq.res_nr.append(deepcopy(indexf+1))
                seq.known.append(False)
                seq.plane_nr.append(1000)
                print "%d%s" % (indexf+1, fastaseq),
            print "\n"
                    
        if FASTAFlag==False:
            for l in range(len(linia_seq)):
                if linia_seq[l]!="\n":
                    iii=ResToNr(linia_seq[l].split()[0]) # check whether seq file contains correct aa names 
                    seq.res_type.append(deepcopy(linia_seq[l].split()[0]))
                    seq.res_nr.append(deepcopy(l+1))
                    seq.known.append(False)
                    seq.plane_nr.append(1000)
        for dim in range(klasy.CPlane.main_dim):
            seq.potentially_present[klasy.CPlane.peak_type.nucl_name[dim]]=True
        for spec in specx:
            for pt in spec.peak_type:
                for dim in range(2):
                    seq.potentially_present[pt.nucl_name[dim]]=True







    return





def ReadPeakList(peaklist_list):
    peak = CPeak
    for i in range(len(peaklist_list)):
        with open(peaklist_list.auto_name, "r") as peaklist_file:
            peaklist_lines = peaklist_file.readlines()
            peaklist = []
            for indexl, line in enumerate(peaklist_lines):
                if indexl > 1 :
                    peaklist_split = line.split()


    return



def ReadInputFile():
    return


def CheckOverlap():
    return


def CalcCCRRate():
    return


def WriteCCRRate():
    return


print ("start")
Experiments = ReadExpSet()