#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mar 15 12:55 2022

@author: Paulina Bartosinska-Marzec
"""
"""  """

""" The script to compering peak list """


from copy import deepcopy
import math
import sys
import os

if len(sys.argv)==1:
    sys.exit("""
    For run script type in command line:
        python3 read_point_list_and_compere [peak list1 path] [peak list2 path] [peak list3 path]...[peak listN path] --dim []
                        
        --dim [] - dimentionality of peak list

    additionaly you can add:
        --comp2list [num] [num]\t- to additionaly compare two of peak list from above
        --out_dir [path]\t- path to direction where output file will be put
        --fl [path]\t- path to file with list of relativ path to peaklists - every set of file must separate by ------
        --intens\t- print peak height, no peak possition\n""")



"""Commantline reading"""



if "--fl" in sys.argv:                                     # number of dimentions
    i = sys.argv.index("--fl")
    FileOfList=[]
    print ("Input:")
    with open(sys.argv[i+1], 'r') as fl:
        p_lines = fl.readlines()
        File_Path = []
        for indexl, line in enumerate(p_lines): 
            if "---" in line:
                FileOfList.append(deepcopy(File_Path))
                File_Path.clear()
            else:
                File_Path.append(deepcopy(line[:-1]))
        if len(File_Path)>0:
            FileOfList.append(deepcopy(File_Path))
    for files in FileOfList:
        print (*files, sep="\n")
else:
    File_Path = []
    for indexa, argument in enumerate(sys.argv):
        if indexa>0:
            if os.path.isfile(argument):
                File_Path.append(deepcopy(sys.argv[indexa]))
    print ("Input Files:",*File_Path, sep="\n")

if "--dim" in sys.argv:                                     # number of dimentions
    i = sys.argv.index("--dim")
    Spectra_dim=int(sys.argv[i+1])

if "--comp2list" in sys.argv:                               # order of two list to additional comparison: 0, 1, 2...
    j = sys.argv.index("--comp2list")
    comparison_of_lists = [int(sys.argv[j+1]), int(sys.argv[j+2])]

if "--out_dir" in sys.argv:                                     # output director
    i = sys.argv.index("--out_dir")
    Output_Dir=sys.argv[i+1]
    if not os.path.exists(Output_Dir):
        os.mkdir(Output_Dir)
else:
    Output_Dir = "./"

if "--intens" in sys.argv:    
    print ("Print only intensity")                                
    FlagIntens = True
else: FlagIntens = False






"""Classes"""


class CPeak:
    def __init__(self):
        self.peak_pos = []            # chemical shifts for all nuclei of peak, length depends on dimentionality
        self.peak_intens = 0            # peak height  
        self.descript = ""              # description of peak which is in first column in Sparky-like peak list
        self.aa_number = 0              # the number of the amino acid to which this peak corresponds



"""Reading functions"""

def read_peaklist(peak_list, s_dim,max_sentence_len):
    FloatFlag = False
    with open(peak_list, 'r') as pl:
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
                    # print ("!!! is digit", item[s_dim])
                    for i in range(1,s_dim+1):
                        p_pos.peak_pos.append(deepcopy(int(item[i])))
                else: 
                    # print ("!!! is NOT digit", item[s_dim])
                    FloatFlag = True
                    for i in range(1,s_dim+1):
                        p_pos.peak_pos.append(deepcopy(float(item[i])))
                if len(item)>s_dim+1:
                    p_pos.peak_intens = float(item[s_dim+1])
                p_list.append(deepcopy(p_pos))
                if max_sentence_len < len(p_pos.descript):
                    max_sentence_len = len(p_pos.descript)
                elif max_sentence_len < len(p_pos.peak_pos):
                    max_sentence_len = len(p_pos.peak_pos)
    return p_list, aa_max_number, max_sentence_len, FloatFlag

def compere_peaklist(peak_list, file_path, aa_max_number, max_sentence_len, FloatFlag, output_dir="./"):
    peak_list_name = []
    compere_summary = []
    for i in file_path:
        item = i.split("/")
        p_name = item[len(item)-1][:-5]
        peak_list_name.append(deepcopy(p_name))
        if max_sentence_len < len(p_name):
            max_sentence_len = len(p_name)
    with open ("{}compere_{}.txt".format(output_dir, peak_list_name[0]), 'w') as outputfile:
        print ("{}compere_{}.txt".format(output_dir, peak_list_name[0]))
        print ("File with comparition of:", *peak_list_name, sep=" ", end="\n\n", file=outputfile)
        print ("AA_num {:^{sentence_len}}".format("Description", sentence_len=max_sentence_len), end=" ", file=outputfile)
        for j in peak_list_name:
            print ('{:^{sentence_len}}'.format(str(j), sentence_len=max_sentence_len), end=" ", file=outputfile)
        try:
            print ("Positions are the same?\t{},{} are the same?".format(comparison_of_lists[0],comparison_of_lists[1]),sep=" ", end="\n\n", file=outputfile)
            # print ("AA_num", "Description", *peak_list_name, "Positions are the same?", "{},{} are the same?".format(comparison_of_lists[0],comparison_of_lists[1]),sep=" ", end="\n\n", file=outputfile)
        except NameError:
            print ("Positions are the same?", sep=" ", end="\n\n", file=outputfile)
        Flag_intens = False
        for aa_num in range(1, aa_max_number):
            # print (aa_num)
            des = "None"
            one_aa = ["None"]*len(peak_list)          # a list of length equal to the number of compared lists
            int_aa = ["-"]*len(peak_list)          # a list of length equal to the number of compared lists
            # print (peak_list)
            
            for indexl, list in enumerate(peak_list):       # loop through peak lists
                for indexp, peak in enumerate(list):
                    if peak.aa_number == aa_num:
                        one_aa[indexl]=peak.peak_pos
                        des = peak.descript
                        if peak.peak_intens:
                            # Flag_intens = True
                            int_aa[indexl]=peak.peak_intens
            # print (one_aa)
            no_peaks = one_aa.count("None")
            half_len_one_aa = math.ceil(len(one_aa)/2)
            # print (no_peaks)
            if no_peaks > half_len_one_aa:
                compere_result = "None"
            elif no_peaks < half_len_one_aa:
                if FloatFlag == True:
                    if all(abs(x[j]-one_aa[0][j])<0.6 for x in one_aa for j in range(len(one_aa))): # it's check difference between all dimentions from any peak with first peak
                        compere_result = "OK"
                    else:
                        compere_result = "Change"
                else:
                    if all(x == one_aa[0] for x in one_aa): # it's check that any peak the same like first peak
                        compere_result = "OK"
                    else:
                        compere_result = "Change"

            print ('{:^6} {:^{sentence_len}}'.format(aa_num, des,sentence_len=max_sentence_len), end=" ", file=outputfile)
            for k in range(len(peak_list)):
                if Flag_intens==True:
                    print ('{:^{sentence_len}} {:^14}'.format(str(one_aa[k]), int_aa[k], sentence_len=max_sentence_len), end=" ", file=outputfile)
                else: 
                    print ('{:^{sentence_len}}'.format(str(one_aa[k]), sentence_len=max_sentence_len), end=" ", file=outputfile)
            try:
                if FloatFlag == True:
                    if one_aa[comparison_of_lists[0]]=="None" or one_aa[comparison_of_lists[0]]=="None":
                        print ("{:^23} {:^17}".format(compere_result," "), end="\n", file=outputfile)
                    elif all(abs(float(one_aa[comparison_of_lists[0]][j])- float(one_aa[comparison_of_lists[1]][j]))<0.6 for j in range(len(one_aa))): # it's check difference between all dimentions from any peak with first peak
                        print ("{:^23} {:^17}".format(compere_result,"yes"), end="\n", file=outputfile)
                    else:
                        print ("{:^23} {:^17}".format(compere_result,"no"), end="\n", file=outputfile)
                else:
                    if one_aa[comparison_of_lists[0]]=="None" or one_aa[comparison_of_lists[0]]=="None":
                        print ("{:^23} {:^17}".format(compere_result," "), end="\n", file=outputfile)
                    elif one_aa[comparison_of_lists[0]]== one_aa[comparison_of_lists[1]]:
                        print ("{:^23} {:^17}".format(compere_result,"yes"), end="\n", file=outputfile)
                    else:
                        print ("{:^23} {:^17}".format(compere_result,"no"), end="\n", file=outputfile)
            except NameError: 
                print ("{:^23}".format(compere_result), end="\n", file=outputfile)
            compere_summary.append(deepcopy(compere_result))
                    


    return compere_summary



def print_all_peaklist(peak_list, file_path, aa_max_number, max_sentence_len, output_dir="./"):
    peak_list_name = []
    compere_summary = []
    for i in file_path:
        item = i.split("/")
        p_name = item[len(item)-1][:-5]
        peak_list_name.append(deepcopy(p_name))
        if max_sentence_len < len(p_name):
            max_sentence_len = len(p_name)
    with open ("{}compere_{}.txt".format(output_dir, peak_list_name[0]), 'w') as outputfile:
        print ("{}compere_{}.txt".format(output_dir, peak_list_name[0]))
        print ("File with comparition of:", *peak_list_name, sep=" ", end="\n\n", file=outputfile)
        print ("AA_num {:^{sentence_len}}".format("Description", sentence_len=max_sentence_len), end=" ", file=outputfile)
        for j in peak_list_name:
            print ('{:^{sentence_len}}'.format(str(j), sentence_len=max_sentence_len), end=" ", file=outputfile)
        # Flag_intens = False
        for aa_num in range(1, aa_max_number):
            # print (aa_num)
            des = "None"
            one_aa = ["None"]*len(peak_list)          # a list of length equal to the number of compared lists
            int_aa = ["-"]*len(peak_list)          # a list of length equal to the number of compared lists
            # print (peak_list)
            
            for indexl, list in enumerate(peak_list):       # loop through peak lists
                for indexp, peak in enumerate(list):
                    if peak.aa_number == aa_num:
                        one_aa[indexl]=peak.peak_pos
                        des = peak.descript
                        if peak.peak_intens:
                            # Flag_intens = True
                            int_aa[indexl]=peak.peak_intens
                            
            no_peaks = one_aa.count("None")
            half_len_one_aa = math.ceil(len(one_aa)/2)

            print ('{:^6} {:^{sentence_len}}'.format(aa_num, des,sentence_len=max_sentence_len), end=" ", file=outputfile)
            for k in range(len(peak_list)):
                print ('{:^14}'.format(int_aa[k]), end=" ", file=outputfile)
            print ("\n", file=outputfile)
                    
    return 


"""                      MAIN PROGRAM                      """
if "FileOfList" in globals():
    print ("\nOutput Files:")
    Compere_Set = []
    for File_Set in FileOfList:
        Peak_Lists = []
        AA_max_number = 1
        Max_sentence_len = 0
        for file_name in File_Set:
            p_list, aa_max, Max_sentence_len, FloatFlag = read_peaklist(file_name,Spectra_dim,Max_sentence_len)
            Peak_Lists.append(deepcopy(p_list))
            if aa_max>AA_max_number:
                AA_max_number=aa_max
        try:
            if FlagIntens == True:
                print_all_peaklist(Peak_Lists,File_Path,AA_max_number,Max_sentence_len)
            else: Compere_Set.append(deepcopy(compere_peaklist(Peak_Lists,File_Set,AA_max_number,Max_sentence_len,FloatFlag,Output_Dir)))
        except:
            item = File_Set[0].split("/")
            Output_Dir = ""
            for i in range(len(item)-1):
                Output_Dir += item[i]+"/"
            if FlagIntens == True:
                print_all_peaklist(Peak_Lists,File_Path,AA_max_number,Max_sentence_len)
            else: Compere_Set.append(deepcopy(compere_peaklist(Peak_Lists,File_Set,AA_max_number,Max_sentence_len,FloatFlag,Output_Dir)))
else:
    Peak_Lists = []
    AA_max_number = 1
    Max_sentence_len = 0
    for file_name in File_Path:
        p_list, aa_max, Max_sentence_len, FloatFlag = read_peaklist(file_name,Spectra_dim,Max_sentence_len)
        Peak_Lists.append(deepcopy(p_list))
        if aa_max>AA_max_number:
            AA_max_number=aa_max
    print ("\nOutput Files:")
    
    if FlagIntens == True:
        print_all_peaklist(Peak_Lists,File_Path,AA_max_number,Max_sentence_len,Output_Dir)
    else:
        Result = compere_peaklist(Peak_Lists,File_Path,AA_max_number,Max_sentence_len,FloatFlag,Output_Dir)
    
