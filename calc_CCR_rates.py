#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

"""
Created on Jan 27 11:18 2022

@author: Paulina Bartosinska-Marzec
"""
""" The script to reading Sparky peak list and calculating CCR rates"""
    

import os
import sys
import csv    
import math
from math import atanh
from copy import deepcopy
from typing import Union
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from scipy import stats
from abc import ABC, abstractmethod

from BioNMRdictionary import Res1to3, Res3to1, CCRname2PrettyRateNamePLT
from CCR_dict import CCR_dict, CCRname2Ratename
import argparse #https://docs.python.org/3/library/argparse.html
import logging
import json

# if len(sys.argv)==1:
#     sys.exit("""
#     The script to calculate CCR rates. In directory has to contains:
#     - peak lists (with peak position in points and peak height)
#     - experiment_set
#     - sequence in FASTA format

#     For run script type in command line:
#         python3 calc_CCR_rate.py [file directory]

#     additionally you can add:
#         --seq  - if name of file with amino acid sequence is not 'seq', add this with name of file
#         --refgamma  - if you have file with reference values of CCR rates, add this with name of file
#                       (file must be .csv, columns name should be: AA, psi_angle (or/and phi_angle), CCR name
#         --expset    - if you want use experiments setup file with different filename than "experiments_set.txt" (structure of file must be the same as original file) """)


"""Command line reading"""


parser = argparse.ArgumentParser(
                    prog='calc_CCR_rates',
                    description="""
The script to calculate cross-correlated relaxation (CCR) rates, using peak intensities from the two spectra, reference and transfer (quantitative gamma approach). 
   
Our scripts can work with any CCR rates, but a few of them are pre-defined (see CCR_dict.py). 

Files required:  
- peak lists with peak positions in ppm and peak heights (`name+"_ppm.list"`)
- peak lists with peak positions in spectral points and peak heights (`name+"_points.list"`)
- JSON file with experiments description (default: `input.json`)   
- sequence in FASTA format (default: `seq`)     

Additional files: 
- peak lists with peak names and uncertainties (`name+"_peaks_noise.list"`)
""",
                    epilog='Text at the bottom of help')

parser.add_argument("file_directory", metavar="file_directory", type=Path, 
                    help="path to directory with all required files")

parser.add_argument("-s", "--seq", dest='seq_file_name', type=Path, 
                    help="name of the file with amino acid sequence (default: `seq`)")

parser.add_argument("-pld", "--peaklist_dir", type=Path, help="directory with all input peak lists (default: identical as file_directory)")

parser.add_argument("-r", "--refgamma", type=Path, 
                    help=argparse.SUPPRESS) #"if you have a file with reference values of CCR rates, add this with the name of the file (file must be .csv, columns names should be: AA, CCR_name_1, CCR_name_2, CCR_name_2, ...)")

parser.add_argument("-e", "--expset", type=Path, 
                    help="experiments setup file (default: `input.json`)")

parser.add_argument("-pub", "--publication", dest='PublicationFlag', action="store_true",  default=False,
                    help=argparse.SUPPRESS) #"if you want the picture outputs to be in the publication size")

parser.add_argument("-pres", "--presentation", dest='PresentationFlag', action="store_true",  default=False,
                    help=argparse.SUPPRESS)#"if you want the picture outputs to be in the presentation size")

args = parser.parse_args()

file_directory = os.path.realpath(args.file_directory)
print ("\nFile directory",file_directory)

if args.seq_file_name:
    seq_file_name = args.seq_file_name 
else: 
    seq_file_name = file_directory+"/seq"
    if not os.path.exists(seq_file_name):
        print (f"there is no file: {seq_file_name}")
print("File with amino acid sequence: {}".format(seq_file_name)) 


if args.refgamma: 
    refgammaFlag = True
    gamma_cal_file_name = args.refgamma
    print ("Reference CCR rates are included from: {}".format(gamma_cal_file_name)) 
else: 
    refgammaFlag = False
    gamma_cal_file_name = None


if args.expset:
    exp_file_name = "{}/{}".format(file_directory,args.expset)
    print("Experiment set description file: {}".format(exp_file_name)) 
else: exp_file_name = "{}/input.json".format(file_directory)

if args.peaklist_dir:
    peaklist_directory = args.peaklist_dir
else:
    peaklist_directory = file_directory
print(f"The peak lists are taken from the directory: {peaklist_directory}")


RaportDir = f"{file_directory}/all_outputs/"
RaportBoxFile = f"{RaportDir}RaportBox.txt"
if not os.path.exists(RaportDir):
        os.mkdir(RaportDir)

if args.PublicationFlag:
    output_style = "publication"
elif args.PresentationFlag:
    output_style = "presentation"
else:
    output_style = "none"


# 

aminoacids=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    


"""Classes"""

class CSeq:
    def __init__(self):
        self.aa_name = ""
        self.aa_number = 0 

class CCRSet:
    
    def __init__(self,file_directory:Path,peaklist_directory:Path, exp_file:Path,seq_dict:dict,ref_flag=False):
        self.__working_dir = file_directory
        self.__peaklist_dir = peaklist_directory
        self.ccr_set = []        # type: list[CCRClass]
        self.__protein_seq = seq_dict
        self.to_compere_dict = {}
        self.__ref_flag = ref_flag
        
        self.ReadExpSet_DICT_JSON(exp_file)
        for exp in self.ccr_set:
            exp.read_peak_files(self.__peaklist_dir,self.__protein_seq)


    def ReadExpSet_DICT_JSON(self,exp_file):               # wczytywanie danych z pliku experiments_set.txt do klasy CCRExp
        with open(exp_file, "r") as json_data_file:
            print(exp_file)
            json_data = json.load(json_data_file)

            for exp_key, exp_properties in json_data.items():
                # exp_properties = json_data[properties]
                # if exp_properties["type_of_CCR"] not in CCR_dict:
                #     exp_name = exp_properties["type_of_CCR"]
                if "symmetrical_reconversion" in exp_properties and exp_properties["symmetrical_reconversion"]==True:
                    one_experiment = CCR_SymRec(exp_properties)
                elif len(exp_properties["NS"])==4:
                    one_experiment = CCR_SymRec(exp_properties)
                else:
                    one_experiment = CCR_normal(exp_properties)

                if exp_properties["type_of_CCR"] not in self.to_compere_dict:
                    self.to_compere_dict[exp_properties["type_of_CCR"]] = [len(self.ccr_set)]
                else:
                    self.to_compere_dict[exp_properties["type_of_CCR"]].append(deepcopy(len(self.ccr_set)))
                
                self.ccr_set.append(deepcopy(one_experiment))
                # print ("{} - {}\n\t reference exp: {}\n\t transfer exp: {}".format(exp_properties["type_of_CCR"],exp_properties["other"],
                #                                                                    one_experiment.ref_name,one_experiment.trans_name))
#                 RaportBox.write(f"""
# Experiment number: {len(self.ccr_set)}
# CCR name: {exp_properties["type_of_CCR"]}
# n dim: {int(exp_properties["dimension"])}
# angle pos: {CCR_dict[exp_properties["type_of_CCR"]]["angle_pos"]}
# ns: {exp_properties["NS"]}
# Tc: {float(exp_properties["TC"]) }
# other: {exp_properties["other"]}

# """)
    def Prepare_peaklist(self, seq_dict:dict,):
        for aminoacids_number in seq_dict:
            res = CResidue(aa_num = aminoacids_number,
                           aa_name = seq_dict[aminoacids_number])
            self._peaks.append(deepcopy(res))

    def CalcCCRrates(self):
        for one_exp in self.ccr_set:
            if one_exp.is_peaklist():
                Add_text = one_exp.Additional_text()
                print_raport(f"\nStarting calculation of CCR rates for: {one_exp.CCRname()}{Add_text}\n")
                one_exp.check_overlap_for_all_peaks()
                one_exp.calc_ccr_rate()
                # print_raport(f"CalcCCRrates: {len(one_exp.peak_list())}")
                # for one_peak in one_exp.peak_list():
                #     print_raport(f"{one_peak.aa_number}: {one_peak.ccr_rate}")
                one_exp.WriteCCRRate_small()
                print_raport("Calculation of CCR rates for {}{} completed".format(one_exp.CCRname(),Add_text))
        self.Write_ALL_CCRRate_CSV()
        

    def Compere_with_reference(self,ref_dict):
        print_raport ("\n=== Comparing with reference values of CCR rates ===\n")
        for one_exp in self.ccr_set:
            if one_exp.is_peaklist():
                one_exp.Add_ref_gamma(ref_dict["seq_num"],ref_dict[one_exp.CCRname()])
                if one_exp.is_reference:
                    print(f"\nComparing CCR rates with calculated for {one_exp.CCRname()}")
                    one_exp.calc_theor_Ix()
                    one_exp.check_if_fatal_error()
                    # one_exp.WriteCCRRate_all_info()
                    # one_exp.WriteCCRRate_all_info_CSV()
                    one_exp.WriteCCRRate_all_info_new(file_type='csv')
                    one_exp.WriteCCRRate_all_info_new(file_type='txt')
                    one_exp.plot_gamma_gamma(style=output_style)
                    one_exp.plot_error_histogram()
        # self.plot_gamma_gamma_all_together(self.ccr_set, style=output_style)
        print_raport ("=== Comparing with reference values of CCR rates completed ===\n")

    def Compere_diff_ver_exp(self):
        if self.__ref_flag:
            if any(len(self.to_compere_dict[CCR_type])>1 for CCR_type in self.to_compere_dict):
                print_raport ("\n=== Comparing diffrent versions of CCR rates with reference values ===\n")
                for CCR_type in self.to_compere_dict:
                    if len(self.to_compere_dict[CCR_type])>1:
                        ExpTable = []
                        for exp_number in self.to_compere_dict[CCR_type]:
                            ExpTable.append(deepcopy(self.ccr_set[exp_number]))
                        self.plot_gamma_gamma_and_diff_together(ExpTable)
                        self.plot_cross_intens_theor_exp_together(ExpTable)
                        self.plot_error_histogram_together(ExpTable)
                        print_raport (f"\n{CCR_type} - completed\n")
                print_raport ("=== Comparing different versions of CCR rates with reference values completed ===\n")
        else:
            if any(len(self.to_compere_dict[CCR_type])>1 for CCR_type in self.to_compere_dict):
                print_raport ("\n=== Comparing experiments with different numbers of NUS points ===\n")
                for CCR_type in self.to_compere_dict:
                    if len(self.to_compere_dict[CCR_type])>1:
                        max_number = 0.0
                        min_number = 10000000.0
                        min_exp = -1
                        max_exp = -1 
                        for exp_number in self.to_compere_dict[CCR_type]:
                            nus_number = self.ccr_set[exp_number].Additional_text()[1:-3]
                            if nus_number.isdigit() == True:
                                nus_number = float(nus_number)
                                if nus_number > max_number:
                                    max_number = nus_number
                                    max_exp = exp_number
                                if nus_number < min_number:
                                    min_number = nus_number
                                    min_exp = exp_number
                        self.plot_gamma_exp_vs_epx(self.ccr_set[min_exp],self.ccr_set[max_exp], str(int(min_number)),str(int(max_number)))
                        print_raport (f"\n{CCR_type} - completed\n")
                print_raport ("\n=== Comparing experiments with different numbers of NUS points completed ===\n")
    
    @staticmethod
    def prepare_gamma_calc_vs_gamma_exp_table(exp_tab:list[CCRClass],
                                              min_max_value:list) -> tuple[list[list[float]],list[Union[int,float]]]:
        set_of_data = []
        
        for indext in range(len(exp_tab)):
            gamma_calculated = []
            gamma_experimental = []
            gamma_calc_error = []
            fatal_error_point_calc = []
            fatal_error_point_exp = []
            fatal_error_point_err = []
            for one_peak in exp_tab[indext].peak_list():
                if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
                    if one_peak.fatal_error:
                        fatal_error_point_calc.append(deepcopy(one_peak.gamma_ref))
                        fatal_error_point_exp.append(deepcopy(one_peak.ccr_rate,))
                        fatal_error_point_err.append(deepcopy(one_peak.ccrrate_error_value))
                    else:
                        gamma_calculated.append(deepcopy(one_peak.gamma_ref))
                        gamma_experimental.append(deepcopy(one_peak.ccr_rate))
                        gamma_calc_error.append(deepcopy(one_peak.ccrrate_error_value))
                    min_max_value = check_if_min_max(one_peak.gamma_ref,one_peak.ccr_rate,min_max_value)
            set_of_data.append(deepcopy([gamma_calculated, gamma_experimental, gamma_calc_error,
                                         [fatal_error_point_calc,fatal_error_point_exp,fatal_error_point_err]]))
        return set_of_data,min_max_value
    
    @staticmethod
    def prepare_gamma_calc_vs_gamma_exp_table_diff_min_max(exp_tab:list[CCRClass]) -> tuple[list[list[float]],list[list[Union[int,float]]]]:
        set_of_data = []
        set_of_min_max_value = []
        for indext in range(len(exp_tab)):
            gamma_calculated = []
            gamma_experimental = []
            gamma_calc_error = []
            fatal_error_point_calc = []
            fatal_error_point_exp = []
            fatal_error_point_err = []
            min_max_value = [+100.0,-100.0] #minimal and maximal value of CCR rate or reference gamma - it is necessary to make plot  
            for one_peak in exp_tab[indext].peak_list():
                if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
                    if one_peak.fatal_error:
                        fatal_error_point_calc.append(deepcopy(one_peak.gamma_ref))
                        fatal_error_point_exp.append(deepcopy(one_peak.ccr_rate))
                        fatal_error_point_err.append(deepcopy(one_peak.ccrrate_error_value))
                    else:
                        gamma_calculated.append(deepcopy(one_peak.gamma_ref))
                        gamma_experimental.append(deepcopy(one_peak.ccr_rate))
                        gamma_calc_error.append(deepcopy(one_peak.ccrrate_error_value))
                    min_max_value = check_if_min_max(one_peak.gamma_ref,one_peak.ccr_rate,min_max_value)
            set_of_data.append(deepcopy([gamma_calculated, gamma_experimental, gamma_calc_error,
                                         [fatal_error_point_calc,fatal_error_point_exp,fatal_error_point_err]]))
            set_of_min_max_value.append(deepcopy(min_max_value))
        return set_of_data,set_of_min_max_value
    
    @staticmethod
    def prepare_gamma_calc_vs_gamma_exp_table_general(exp_tab:list[CCRClass],min_max_value=[]) -> tuple[list[dict[list[float]]],list[list[Union[int,float]]]|list[float,float]]:
        if len(min_max_value) == 0:
            minmaxFlag = True
        else:
            minmaxFlag = False
        set_of_data = []
        set_of_min_max_value = []
        for indext in range(len(exp_tab)):
            gamma_calculated = []
            gamma_experimental = []
            gamma_calc_error = []
            fatal_error_point_calc = []
            fatal_error_point_exp = []
            fatal_error_point_err = []
            if minmaxFlag:
                min_max_value = [+100.0,-100.0] #minimal and maximal value of CCR rate or reference gamma - it is necessary to make plot  
            for one_peak in exp_tab[indext].peak_list():
                if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
                    if one_peak.fatal_error:
                        fatal_error_point_calc.append(deepcopy(one_peak.gamma_ref))
                        fatal_error_point_exp.append(deepcopy(one_peak.ccr_rate))
                        fatal_error_point_err.append(deepcopy(one_peak.ccrrate_error_value))
                    else:
                        gamma_calculated.append(deepcopy(one_peak.gamma_ref))
                        gamma_experimental.append(deepcopy(one_peak.ccr_rate))
                        gamma_calc_error.append(deepcopy(one_peak.ccrrate_error_value))
                    min_max_value = check_if_min_max(one_peak.gamma_ref,one_peak.ccr_rate,min_max_value)
            set_of_data.append(deepcopy({"good":[gamma_calculated, gamma_experimental, gamma_calc_error],
                                         "fatal":[fatal_error_point_calc,fatal_error_point_exp,fatal_error_point_err]}))
            if minmaxFlag:
                set_of_min_max_value.append(deepcopy(min_max_value))
        
        if minmaxFlag:
            return set_of_data,set_of_min_max_value
        else:
            return set_of_data, min_max_value

    @staticmethod
    def prepare_gamma_exp_vs_diff_nus_table(exp_tab:list[CCRClass],
                                            min_max_value:list) -> tuple[list[list[float]],list[str],list[Union[int,float]]]:
        set_of_data = [[],[]]
        
        max_number = 0
        min_number = 10000000
        min_exp = -1
        max_exp = -1 
        for indext, one_exp in enumerate(exp_tab):
            nus_number = one_exp.other_info()[1:-3]
            if nus_number.isdigit() == True:
                nus_number = int(nus_number)
                if nus_number > max_number:
                    max_number = nus_number
                    max_exp = indext
                if nus_number < min_number:
                    min_number = nus_number
                    min_exp = indext
        min_max_discript = [exp_tab[min_exp].other_info(),exp_tab[max_exp].other_info()]
        peaks_table_min = exp_tab[min_exp].peak_list()
        peaks_table_max = exp_tab[max_exp].peak_list()
        for indexp, one_peak in enumerate(peaks_table_max):
            if one_peak.is_ccr_rate and peaks_table_min[indexp].is_ccr_rate and one_peak.aa_name!="G":
                set_of_data[0].append(deepcopy(peaks_table_min[indexp].ccr_rate))
                set_of_data[1].append(deepcopy(one_peak.ccr_rate))
                min_max_value = check_if_min_max(peaks_table_min[indexp].ccr_rate,one_peak.ccr_rate,min_max_value)
            
        return set_of_data,min_max_discript,min_max_value

    def plot_gamma_gamma_all_together(self,exp_tab:list[CCRClass],transparent_plot=False, nrows=3, ncols=3, add2="", style="none"):
        if nrows==3 and ncols==3:
            nrows,ncols = set_nrow_ncol(len(exp_tab))

        colours_set = ['red','gold','darkturquoise','darkviolet','darkorange','limegreen','blue','lightpink']
        fig, axs = plt.subplots(nrows=nrows,ncols=ncols, layout='constrained')
        # if not style:
        #     fig.suptitle(f"Comparision of experimental CCR rates with structure-predicted CCR rates based on:\n\"{gamma_cal_file_name}\" ")
        
        set_of_data,min_max_value = self.prepare_gamma_calc_vs_gamma_exp_table_general(exp_tab)
        # print(set_of_data)
        
        # first row - plot each NUS number in particular plot
        for indexa, ax in enumerate(axs.flat):
            if indexa < len(exp_tab):
                exp = exp_tab[indexa]
                Add_text = exp.Additional_text()
                if exp.is_peaklist() and exp.is_reference():
                    exp_ptl_data = set_of_data[indexa]
                    ax.axline([0,0],slope=1, linestyle=(0, (3, 3)), linewidth=1, color='darkgray', label='y=x') 
                    ax.scatter(exp_ptl_data['good'][0],exp_ptl_data['good'][1],s=2, color='#252525ff', ) #
                    ax.errorbar(exp_ptl_data['good'][0], exp_ptl_data['good'][1], 
                            yerr=exp_ptl_data['good'][2], ecolor='#c2a0f8ff', 
                            fmt='o', markerfacecolor='#252525ff', markeredgecolor='none',markersize=2) #252525ff
                    
                    
                    
                    #fatal errors points:
                    if style != 'none':
                        dot_colours1 = '#252525ff'
                        dot_colours2 = '#c2a0f8ff'
                    else:
                        dot_colours1 = 'red'
                        dot_colours2 = 'red'
                    if len(exp_ptl_data['fatal'][0])>0 and len(exp_ptl_data['fatal'][1])>0 and len(exp_ptl_data['fatal'][2])>0:
                        ax.errorbar(exp_ptl_data['fatal'][0], exp_ptl_data['fatal'][1], 
                                yerr=exp_ptl_data['fatal'][2], ecolor=dot_colours2, 
                            fmt='o', markerfacecolor=dot_colours1, markeredgecolor='none',markersize=2) #
                        ax.scatter(exp_ptl_data['fatal'][0],exp_ptl_data['fatal'][1],s=2, color=dot_colours1) #
                        
                        
                    weighted_reg_dict = WeightedLRegression_expresion_by_hand(x=exp_ptl_data['good'][0],
                                                                            y=exp_ptl_data['good'][1],
                                                                            uncertainty_val=exp_ptl_data['good'][2])
                    
                    label_text = "{}:\t{} peaks,\t{},  delta_a: {:.3f}, delta_b: {:.3f}    r2 = {:.2f}    factor \'a\' = {:.1f}, factor \'b\' = {:.1f}".format(Add_text[1:],
                                                                                                                                                                weighted_reg_dict["equation"], 
                                                                                                                                                                weighted_reg_dict["slope_uncertainty"], 
                                                                                                                                                                weighted_reg_dict["intercept_uncertainty"],
                                                                                                                                                                weighted_reg_dict["r2"],
                                                                                                                                                                weighted_reg_dict["factor_a"],
                                                                                                                                                                weighted_reg_dict["factor_b"],
                                                                                                                                                                len(exp_ptl_data['good'][1]))
                    
                    label_text2 = "{}\ndelta_a: {:.3f}\ndelta_b: {:.3f}\nfactor \'a\' = {:.1f}".format(weighted_reg_dict["equation"], 
                                                                                                                weighted_reg_dict["slope_uncertainty"], 
                                                                                                                weighted_reg_dict["intercept_uncertainty"],
                                                                                                                weighted_reg_dict["factor_a"],)
                    
                    text_pos_x = min_max_value[indexa][0] - (abs(min_max_value[indexa][0])+abs(min_max_value[indexa][1]))/5
                    text_pos_y = min_max_value[indexa][1] + (abs(min_max_value[indexa][0])+abs(min_max_value[indexa][1]))/5
                    
                    if style == "publication":
                        ax.set_title(f'{CCRname2PrettyRateNamePLT(exp.CCRname())}', fontsize=8.5, color ='#252525ff' )
                        trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
                        ax.text(-0.05, 1.0, f'Exp. {exp.CCRname()[4:]}', transform=ax.transAxes + trans,
                                fontsize=7, verticalalignment='top')
                        
                        # ax.text(text_pos_x,text_pos_y,
                        #     f'Exp. {exp.CCRname()[4:]}', 
                        #     horizontalalignment='left', verticalalignment='top',
                        #     fontsize=7)
                        # ax.text(text_pos_x,text_pos_y,
                        #     f'$R^{2}$ = {weighted_reg_dict["r2"]:.2f}', 
                        #     horizontalalignment='left', verticalalignment='top',
                        #     fontsize=7)
                    elif style == "presentation":
                        ax.set_title(f'{CCRname2PrettyRateNamePLT(exp.CCRname())}', fontsize=12, color = colours_set[indexa])
                        # ax.text(text_pos_x,text_pos_y,
                        #     f'$R^{2}$ = {weighted_reg_dict["r2"]:.2f}', 
                        #     horizontalalignment='left', verticalalignment='top',
                        #     fontsize=9)
                    else:
                        ax.axline([0,weighted_reg_dict["intercept"]],slope=weighted_reg_dict["slope"], 
                              linestyle=(0, (5, 5)), linewidth=1, color='mediumorchid', 
                              label=label_text)
                        # ax.text(text_pos_x,text_pos_y,
                        #     f'$R^{2}$ = {weighted_reg_dict["r2"]:.2f}', 
                        #     horizontalalignment='left', verticalalignment='top',
                        #     fontsize=9)
                        ax.text(text_pos_x,text_pos_y,
                            label_text2, 
                            horizontalalignment='left', verticalalignment='top',
                            fontsize=9)
                        ax.set_title(f"{exp.CCRname()} {Add_text[1:]}", fontsize=10) 

                    ax = setup_plot_area(ax,min_max_value[indexa])
                    ax.tick_params(labelsize=8)

        #mediumpurple
        

        if add2 == '':
            add2 = exp_tab[0].Additional_text()

        if style == "publication":
            # plt.subplots_adjust(hspace=0.3,wspace = 0.3)
            fig.supxlabel('structure-predicted \u0393, $s^{-1}$', fontsize=10)
            fig.supylabel('experimental \u0393, $s^{-1}$', fontsize=10)
            # fig.supxlabel('2D conventional \u0393, $s^{-1}$', fontsize=10)
            # fig.supylabel('4D NUS \u0393, $s^{-1}$', fontsize=10)
            cm = 1/2.54  # centimeters in inches
            fig.set_figwidth(17.4*cm)
            if nrows == 3:
                fig.set_figheight(17*cm)
            elif nrows == 2:
                fig.set_figheight(10*cm)
            fig.savefig("{}/{}_all_exp_vs_calc{}.eps".format(file_directory,str(gamma_cal_file_name)[:-4],add2), bbox_inches="tight", 
                   transparent=transparent_plot, dpi=1200, format='eps')  # pad_inches=0.3, 
        elif style == "presentation":
            # plt.subplots_adjust(hspace=0.3,wspace = 0.25)
            fig.supxlabel('structure-predicted \u0393, $s^{-1}$', fontsize=11)
            fig.supylabel('experimental \u0393, $s^{-1}$', fontsize=11)
            lines, labels = axs.flat[0].get_legend_handles_labels()
            fig.legend(lines, labels, loc='lower left',fontsize=10)
            fig.set_figwidth(12)
            fig.set_figheight(6)
            fig.savefig("{}/{}_all_exp_vs_calc{}.svg".format(file_directory,str(gamma_cal_file_name)[:-4],add2), bbox_inches="tight", 
                     transparent=transparent_plot, dpi=1200, format='svg')
        else:
            # plt.subplots_adjust(hspace=0.3,wspace = 0.25)
            fig.supxlabel('structure-predicted \u0393, $s^{-1}$')
            fig.supylabel('experimental \u0393, $s^{-1}$')
            fig.set_figwidth(ncols*4)
            fig.set_figheight(nrows*4)
            # fig.legend()
        fig.savefig("{}/{}_all_exp_vs_calc{}.png".format(file_directory,str(gamma_cal_file_name)[:-4],add2), bbox_inches="tight", 
                transparent=transparent_plot, dpi=1200, format='png')
        plt.close()
        plt.clf()

    def plot_gamma_gamma_and_diff_together(self,exp_tab:list[CCRClass],transparent_plot=False,style='none'):
        ccr_name = exp_tab[0].CCRname()
        # plt.rcParams['font.size'] = '14'
        min_max_value = [+100.0,-100.0]
        fig, axs = plt.subplots(nrows=2,ncols=len(exp_tab),layout='constrained')
        if not style:
            fig.suptitle(f"Comparision of CCR rates ({ccr_name}) for diffrent experiment variant with structure-predicted CCR rates based on:\n\"{gamma_cal_file_name}\"")
        
        set_of_data,min_max_value = self.prepare_gamma_calc_vs_gamma_exp_table_general(exp_tab,min_max_value=min_max_value)
        set_of_data_min_max_NUS,min_max_discript,min_max_value = self.prepare_gamma_exp_vs_diff_nus_table(exp_tab,min_max_value)
    
        # first row - plot each NUS number in particular plot
        for indext, exp in enumerate(exp_tab):
            exp_ptl_data = set_of_data[indext]
            print (exp_ptl_data)
            Add_text = exp.Additional_text()
            axs[0,indext].axline([0,0],slope=1, linestyle=(0, (3, 3)), linewidth=1, color='darkgray') 
            axs[0,indext].scatter(exp_ptl_data['good'][0],exp_ptl_data['good'][1],s=10,) #
            axs[0,indext].errorbar(exp_ptl_data['good'][0], exp_ptl_data['good'][1], 
                        yerr=exp_ptl_data['good'][2], 
                        fmt='none', color='#252525ff') #
            #fatal errors points:
            if len(exp_ptl_data['fatal'][0])>0 and len(exp_ptl_data['fatal'][1])>0 and len(exp_ptl_data['fatal'][2])>0:
                axs[0,indext].scatter(exp_ptl_data['fatal'][-1][0],exp_ptl_data['fatal'][-1][1],s=3, color='red') #
                axs[0,indext].errorbar(exp_ptl_data['fatal'][-1][0], exp_ptl_data['fatal'][-1][1], 
                        yerr=exp_ptl_data['fatal'][-1][2], 
                        fmt='none', color='red') #
            if len(exp_ptl_data['good'][0])>0:
                weighted_reg_dict = WeightedLRegression_expresion_by_hand(exp_ptl_data['good'][0], exp_ptl_data['good'][1], exp_ptl_data['good'][2])
                label_text = "{}:\t{}\ndelta_a: {:.3f}, delta_b: {:.3f}\tr2 = {:.2f}\nfactor \'a\' = {:.1f}, factor \'b\' = {:.1f}".format(Add_text[1:],
                                                                                                                                                                weighted_reg_dict["equation"], 
                                                                                                                                                                weighted_reg_dict["slope_uncertainty"], 
                                                                                                                                                                weighted_reg_dict["intercept_uncertainty"],
                                                                                                                                                                weighted_reg_dict["r2"],
                                                                                                                                                                weighted_reg_dict["factor_a"],
                                                                                                                                                                weighted_reg_dict["factor_b"],
                                                                                                                                                                len(exp_ptl_data['good'][1]))
                axs[0,indext].axline([0,weighted_reg_dict["intercept"]],slope=weighted_reg_dict["slope"], 
                            linestyle=(0, (5, 5)), linewidth=1.5, color='purple', 
                            label=label_text)
                # axs[0,indext].legend()
            axs[0,indext].set_title(Add_text[1:], fontsize=10)
            axs[0,indext] = setup_plot_area(axs[0,indext],min_max_value)
            axs[0,indext].tick_params(labelsize=8)

        # second row
        axs[1,0].axline([0,0],slope=1, linestyle=(0, (3, 3)), linewidth=1, color='darkgray', label='x=y') 
        for indext, exp in enumerate(exp_tab):
            Add_text = exp.Additional_text()
            axs[1,0].scatter(exp_ptl_data['good'][0],
                             exp_ptl_data['good'][1],
                             s=10,label=str(Add_text[1:])) #
        axs[1,0] = setup_plot_area(axs[1,0],min_max_value)
        axs[1,0].tick_params(labelsize=8)
        fig.legend(fontsize="medium",loc='lower left', bbox_to_anchor=(0.3, 0.1))

        for indext in range(1,len(exp_tab)-1):
            axs[1][indext].set_axis_off()

        axs[1,len(exp_tab)-1].axline([0,0],slope=1, linestyle=(0, (3, 3)), linewidth=1, color='darkgray') 
        axs[1,len(exp_tab)-1].scatter(set_of_data_min_max_NUS[1],set_of_data_min_max_NUS[0],s=10,) #
        axs[1,len(exp_tab)-1]= setup_plot_area(axs[1,len(exp_tab)-1],min_max_value)
        axs[1,len(exp_tab)-1].set_xlabel(min_max_discript[1]+' NUS \u0393, $s^{-1}$',fontsize=10)
        axs[1,len(exp_tab)-1].set_ylabel(min_max_discript[0]+' NUS \u0393, $s^{-1}$',fontsize=10)
        axs[1,len(exp_tab)-1].tick_params(labelsize=8)

        fig.supxlabel('structure-predicted \u0393, $s^{-1}$')
        fig.supylabel('experimental \u0393, $s^{-1}$')

        # plt.subplots_adjust(hspace=0.3,wspace = 0.3)

        if style == "publication":
            cm = 1/2.54  # centimeters in inches
            fig.set_figwidth(17.4*cm)
            fig.set_figheight(17*cm)
            fig.savefig("{}/{}_exp_vs_calc_diff.eps".format(file_directory, ccr_name), bbox_inches="tight", 
                     transparent=transparent_plot, dpi=1200, format='eps')
        elif style == "presentation":
            fig.set_figwidth(len(exp_tab)*3.5)
            fig.set_figheight(7)
            fig.savefig("{}/{}_exp_vs_calc_diff.svg".format(file_directory, ccr_name), bbox_inches="tight", 
                     transparent=transparent_plot, dpi=1200, format='svg')
        else:
            fig.set_figwidth(len(exp_tab)*3.5)
            fig.set_figheight(7)
        fig.savefig("{}/{}_exp_vs_calc_diff.png".format(file_directory, ccr_name), bbox_inches="tight", 
                 transparent=transparent_plot, dpi=1200, format='png')   #pad_inches=0.3,
        plt.close()
        plt.clf()

    def plot_gamma_diff(self,exp_tab:list[CCRClass],transparent_plot=False):
        ccr_name = exp_tab[0].CCRname()
        plt.rcParams['font.size'] = '14'
        # Plotting both the curves simultaneously
        plt.axline([0,0],slope=1, linestyle=(0, (5, 5)), linewidth=1.5, color='darkgray', label='x=y')

        min_max_value = [+100.0,-100.0]         #minimal and maximal value of CCR rate or reference gamma - it is necessary to make plot  
        for indexe, one_exp in enumerate(exp_tab):
            Add_text = one_exp.Additional_text()
            gamma_calculated = []
            gamma_experimental = []
            for one_peak in one_exp.peak_list():
                if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
                    gamma_calculated.append(deepcopy(one_peak.gamma_ref))
                    gamma_experimental.append(deepcopy(one_peak.ccr_rate))
                    if one_peak.gamma_ref-one_peak.ccr_rate>5 and indexe==0:
                        plt.annotate(one_peak.aa_number,(one_peak.gamma_ref,one_peak.ccr_rate),
                                    textcoords="offset points",xytext=(0,-12),ha='center')
                    if one_peak.ccr_rate-one_peak.gamma_ref>5 and indexe==0:
                        plt.annotate(one_peak.aa_number,(one_peak.gamma_ref,one_peak.ccr_rate),
                                    textcoords="offset points",xytext=(0,10),ha='center')
                    min_max_value = check_if_min_max(one_peak.gamma_ref,one_peak.ccr_rate,min_max_value)
            plt.scatter(gamma_calculated, gamma_experimental, s=10, label=str(Add_text[indexe][1:]))

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
        plt.savefig("{}/{}_exp_vs_calc_diff.png".format(file_directory, ccr_name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
        # plt.show()
        plt.clf()
        plt.close()

    def plot_cross_intens_theor_exp_together(self,exp_tab:list[CCRClass],transparent_plot=False):
        ccr_name = exp_tab[0].CCRname()
        # plt.rcParams['font.size'] = '14'

        min_max_value = [+100.0,-100.0]
        fig, axs = plt.subplots(ncols=len((exp_tab)))
        fig.suptitle("Comparision of intensity of peaks in transfer spectra {}\n(structure-predicted vs experimental) \nwith diffrent number of NUS points".format(ccr_name))
        
        set_of_data = []
        for indext, one_exp in enumerate(exp_tab):
            Add_text = one_exp.Additional_text()
            Ix_theory = []
            Ix_experimental = []
            for one_peak in exp_tab[indext].peak_list():
                if one_peak.is_peak and one_peak.Ix_theor!=0.0:
                    # print ("there is peak {} for spectra {} ({})".format(one_peak.aa_number,ccr_name,add[indexa]))
                    # print("Intensity of Ix_theory: {}, Ix_experimental: {}, min_max_value: {}".format(one_peak.Ix_theor, one_peak.peak_intens[1],min_max_value))
                    Ix_theory.append(deepcopy(one_peak.Ix_theor))
                    Ix_experimental.append(deepcopy(one_peak.peak_intens[1]))
                    min_max_value = check_if_min_max(one_peak.Ix_theor,one_peak.peak_intens[1],min_max_value)
            set_of_data.append(deepcopy([Ix_theory,Ix_experimental]))
            # print("length of Ix_theory: {}, Ix_experimental: {} for spectra {} ({})".format(len(Ix_theory), len(Ix_experimental),ccr_name,add[indext]))
            # print(Ix_theory, Ix_experimental)

        for indext, one_exp in enumerate(exp_tab):
            Add_text = one_exp.Additional_text()
            axs[indext].axline([0,0],slope=1, linestyle=(0, (3, 3)), linewidth=1, color='darkgray', label='x=y') 
            axs[indext].scatter(set_of_data[indext][0],set_of_data[indext][1],s=10,) #
            axs[indext].set_title(Add_text[indext][1:],fontsize=14)
            axs[indext] = setup_plot_area(axs[indext],min_max_value)
            axs[indext].tick_params(labelsize=8)
            axs[indext].set_aspect('equal', adjustable='box')
        fig.supxlabel('structure-predicted intensity')
        fig.supylabel('experimental intensity')
        fig.set_figwidth(len(exp_tab)*3)
        fig.set_figheight(4)
        # plt.subplots_adjust(hspace=0.3,wspace = 0.3)
        fig.savefig("{}/{}_exp_vs_calc_Ix.png".format(file_directory, ccr_name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
        plt.close()
        plt.clf()

    def plot_error_histogram_together(self,exp_tab:list[CCRClass],transparent_plot=False): # ERROR GAMMA INFO, 9.11.2023
        ccr_name = exp_tab[0].CCRname()
        set_of_data = []
        for indext, one_exp in enumerate(exp_tab):
            Add_text = one_exp.Additional_text()
            gamma_experimental = []
            gamma_calc_error = []            
            peak_intens_auto = []         
            peak_intens_cross = []
            for one_peak in one_exp.peak_list():
                if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G" and one_peak.is_peak:
                    gamma_experimental.append(deepcopy(one_peak.ccr_rate))
                    gamma_calc_error.append(deepcopy(one_peak.ccrrate_error_value)) 
                    peak_intens_auto.append(deepcopy(one_peak.peak_intens[0]))
                    peak_intens_cross.append(deepcopy(one_peak.peak_intens[1]))
            set_of_data.append(deepcopy({"gamma_experimental":gamma_experimental,
                                        "gamma_calc_error":gamma_calc_error,
                                        "peak_intens_auto":peak_intens_auto,
                                        "peak_intens_cross":peak_intens_cross}))

        fig, axs = plt.subplots(nrows=3,ncols=len((exp_tab)))    
        plt.subplots_adjust(hspace=0.4)
        fig.set_figheight(15)
        fig.set_figwidth(6*len((exp_tab)))
        
        for indext, one_exp in enumerate(exp_tab):
            Add_text = one_exp.Additional_text()
            axs[0,indext].scatter(set_of_data[indext]["gamma_experimental"], set_of_data[indext]["gamma_calc_error"], s=5, color='#0066ffff')
            axs[0,indext].set_title("{} {}\n Experimental \u0393 vs Uncertainty value".format(ccr_name,str(Add_text[indext][1:])))
            axs[0,indext].set_ylabel('Uncertainty value')
            axs[0,indext].set_xlabel('Experimental \u0393, $s^{-1}$')

            axs[1,indext].scatter(set_of_data[indext]["peak_intens_auto"], set_of_data[indext]["gamma_calc_error"], s=5, color='#0066ffff')
            axs[1,indext].set_title("Peak intens in reference vs Uncertainty value")
            axs[1,indext].set_ylabel('Uncertainty value')
            axs[1,indext].set_xlabel('Peak intens in reference')

            axs[2,indext].scatter(set_of_data[indext]["peak_intens_cross"], set_of_data[indext]["gamma_calc_error"], s=5, color='#0066ffff')
            axs[2,indext].set_title("Peak intens in transfer vs Uncertainty value")
            axs[2,indext].set_ylabel('Uncertainty value')
            axs[2,indext].set_xlabel('Peak intens in transfer')

        plt.savefig("{}Hist_{}_exp_vs_Uncertainty_all.png".format(RaportDir, ccr_name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
        # plt.show()
        plt.clf()
        plt.close()

    def plot_gamma_exp_vs_epx(self,exp_min,exp_max,min_add,max_add,transparent_plot=False):
        ccr_name = exp_min.CCRname()
        peaks_min = exp_min.peak_list()
        peaks_max = exp_max.peak_list()

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
        plt.savefig("{}/{}_minNUS_vs_maxNUS.png".format(file_directory, ccr_name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
        plt.clf()
        plt.close()

    def Write_ALL_CCRRate_CSV(self):
        new_list = f"{file_directory}/CCRrates.csv"
        with open(new_list, mode='w', newline='') as csv_file:
            headers = ['AA']
            for one_exp in self.ccr_set:
                Add_textt=one_exp.Additional_text()
                headers.append(deepcopy('{}{}'.format(one_exp.CCRname(),Add_textt)))
                headers.append(deepcopy('{}_error'.format(one_exp.CCRname())))            # ERROR GAMMA INFO, 9.11.2023
            writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
            writer.writeheader()
            zero_row = {}
            zero_row['AA'] = 'res'
            for one_exp in self.ccr_set:
                Add_textt=one_exp.Additional_text()
                zero_row['{}{}'.format(one_exp.CCRname(),Add_textt)] = '{}'.format(CCRname2Ratename(one_exp.CCRname()))
            writer.writerow(zero_row)
            for seq_number in self.__protein_seq:
                one_row = {}
                one_row["AA"] = str(seq_number)+Res1to3(self.__protein_seq[seq_number])
                for one_exp in self.ccr_set:
                    Add_textt=one_exp.Additional_text()
                    for one_peak in one_exp.peak_list():
                        if one_peak.aa_number == seq_number:
                            if one_peak.is_ccr_rate == True and one_peak.ccrrate_calculation_error == False and one_peak.aa_name != "G":
                                one_row['{}{}'.format(one_exp.CCRname(),Add_textt)] = '{:.4f}'.format(one_peak.ccr_rate)
                                one_row['{}_error'.format(one_exp.CCRname())] = '{:.4f}'.format(one_peak.ccrrate_error_value)            # ERROR GAMMA INFO, 9.11.2023
                            else:
                                one_row['{}{}'.format(one_exp.CCRname(),Add_textt)] = "nan"
                                one_row['{}_error'.format(one_exp.CCRname())] = "nan"            # ERROR GAMMA INFO, 9.11.2023
                writer.writerow(one_row)
        print ("All CCR rates are in:", new_list)
        return

    def Write_zeroCCRrates_Analisis_CSV(self):
        new_list = f"{RaportDir}zeroCCRrates_Analisis.csv"
        
        s_dim = self.ccr_set[0].s_dim()
        error_residue_dict={}
        error_experiment_dict={}
        error_aa_dict = {}
        invisible_cross_peaks = {}
        aa_dict_nun = {}

        with open(new_list, mode='w', newline='') as csv_file:
            headers = ['CCR name','AA']
            for i in range(1,s_dim+1):
                headers.append(deepcopy('w{}'.format(i)))
            headers.extend(deepcopy(['intensity in reference','Noise in reference','intensity in transfer','Noise in transfer','CCR rate','Uncertainty','Comments','Reference Gamma', 'Ix Theor', 'Ix/Ia', 'Ix/Ia without NS']))
            writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
            writer.writeheader()
            for one_exp in self.ccr_set:
                CCR_name =  one_exp.CCRname()+one_exp.Additional_text()
                if one_exp.is_peaklist() == True:
                    writer.writerow({})
                    invisible_cross_peaks[CCR_name]=[0,0]
                    error_experiment_dict[CCR_name]=0
                    for peak_number, one_peak in enumerate(one_exp.peak_list()):
                        if one_peak.aa_name in aa_dict_nun:
                            aa_dict_nun[one_peak.aa_name]+=1
                        else:
                            aa_dict_nun[one_peak.aa_name]=1                    
                        one_exp.deepVisibilityAnalisis(one_peak,invisible_cross_peaks)

                        if one_peak.is_ccr_rate and one_peak.aa_name != "G" and one_peak.is_gamma_calc:
                            if (one_peak.gamma_ref - one_peak.ccr_rate)/one_peak.ccr_rate >=2 or (one_peak.gamma_ref - one_peak.ccr_rate)/one_peak.ccr_rate <= -2:
                                one_exp.deepResidueAnalisys(writer,peak_number,
                                                            error_residue_dict,
                                                            error_experiment_dict,
                                                            error_aa_dict)

            writer_row = csv.writer(csv_file,delimiter=",")
            writer_row.writerow("")
            
            self.deepExpAnalisys(writer_row,
                    error_residue_dict,
                    error_experiment_dict,
                    error_aa_dict,
                    invisible_cross_peaks,
                    aa_dict_nun)
            
        # print ("All info from peaks for {} rates are in file: {}".format(ccr_name, new_list), file=RaportBox)
        RaportBox.write("Information about very small CCR rates are in file: {}\n".format(new_list))
        return

    def deepExpAnalisys(self,
                        writer_to_file_by_dict:csv,
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
        writer_to_file_by_dict.writerow(["Number of peaks in transfer where:","CCR zero error","intesity < 10*noise level","intesity < noise level"])
        for one_exp in self.ccr_set:
            CCR_name =  one_exp.CCRname()+one_exp.Additional_text()
            if CCR_name in error_experiment_dict and CCR_name in invisible_cross_peaks:
                writer_to_file_by_dict.writerow([CCR_name,error_experiment_dict[CCR_name],invisible_cross_peaks[CCR_name][0],invisible_cross_peaks[CCR_name][1]])
        writer_to_file_by_dict.writerow("")
        writer_to_file_by_dict.writerow(["AA name","CCR zero error"])
        for aa in aminoacids:
            if Res3to1(aa) in error_aa_dict:
                # print_raport(f"{aa},{error_aa_dict[Res3to1(aa)]},{aa_dict_nun[Res3to1(aa)]}\n")
                percent = error_aa_dict[Res3to1(aa)]/aa_dict_nun[Res3to1(aa)]*100
                writer_to_file_by_dict.writerow([aa,error_aa_dict[Res3to1(aa)],f"{percent:.1f} %"])

    def Write_ErrorCCRrates_Analisis_CSV(self):
        new_list = f"{RaportDir}ErrorCCRrates_Analisis.csv"
        
        s_dim = self.ccr_set[0].s_dim()
        error_residue_dict={}
        error_experiment_dict={}
        error_aa_dict = {}
        invisible_cross_peaks = {}
        aa_dict_nun = {}

        with open(new_list, mode='w', newline='') as csv_file:
            headers = ['CCR name','AA']
            for i in range(1,s_dim+1):
                headers.append(deepcopy('w{}'.format(i)))
            headers.extend(deepcopy(['intensity in reference','Noise in reference','intensity in transfer','Noise in transfer','CCR rate','Uncertainty','Comments','Reference Gamma', 'Ix Theor', 'Ix/Ia', 'Ix/Ia without NS']))
            writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
            writer.writeheader()
            for one_exp in self.ccr_set:
                CCR_name =  one_exp.CCRname()+one_exp.Additional_text()
                if one_exp.is_peaklist() == True:
                    writer.writerow({})
                    invisible_cross_peaks[CCR_name] = [0,0]
                    error_experiment_dict[CCR_name] = 0
                    for peak_number, one_peak in enumerate(one_exp.peak_list()):
                        if one_peak.aa_name in aa_dict_nun:
                            aa_dict_nun[one_peak.aa_name]+=1
                        else:
                            aa_dict_nun[one_peak.aa_name]=1
                        one_exp.deepVisibilityAnalisis(one_peak,invisible_cross_peaks)

                        if one_peak.is_ccr_rate and one_peak.aa_name != "G" and one_peak.gamma_ref > one_peak.ccr_rate*5 and one_peak.is_gamma_calc:
                                one_exp.deepResidueAnalisys(writer,peak_number,
                                                        error_residue_dict,
                                                        error_experiment_dict,
                                                        error_aa_dict)

            writer_row = csv.writer(csv_file,delimiter=",")
            writer_row.writerow("")
            
            self.deepExpAnalisys(writer_row,
                    error_residue_dict,
                    error_experiment_dict,
                    error_aa_dict,
                    invisible_cross_peaks,
                    aa_dict_nun)
            
            
        # print ("All info from peaks for {} rates are in file: {}".format(ccr_name, new_list), file=RaportBox)
        RaportBox.write("Information about CCR rates which are off are in file: {}\n".format(new_list))
        return

    def Write_diffCCRrates_Analisis_CSV(self):
        new_list = f"{RaportDir}diffCCRrates_Analisis.csv"
        
        s_dim = self.ccr_set[0].s_dim()
        error_residue_dict={}
        error_experiment_dict={}
        error_aa_dict = {}
        invisible_cross_peaks = {}
        aa_dict_nun = {}

        with open(new_list, mode='w', newline='') as csv_file:
            headers = ['CCR name','AA']
            for i in range(1,s_dim+1):
                headers.append(deepcopy('w{}'.format(i)))
            headers.extend(deepcopy(['intensity in reference','Noise in reference','intensity in transfer','Noise in transfer','CCR rate','Uncertainty','Comments','Reference Gamma', 'Ix Theor', 'Ix/Ia', 'Ix/Ia without NS']))
            writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
            writer.writeheader()
            for one_exp in self.ccr_set:
                CCR_name =  one_exp.CCRname()+one_exp.Additional_text()
                if one_exp.is_peaklist() == True:
                    writer.writerow({})
                    invisible_cross_peaks[CCR_name]=[0,0]
                    error_experiment_dict[CCR_name]=0
                    for peak_number, one_peak in enumerate(one_exp.peak_list()):
                        if one_peak.aa_name in aa_dict_nun:
                            aa_dict_nun[one_peak.aa_name]+=1
                        else:
                            aa_dict_nun[one_peak.aa_name]=1
                        one_exp.deepVisibilityAnalisis(one_peak,invisible_cross_peaks)

                        if one_peak.is_ccr_rate and one_peak.aa_name != "G" and one_peak.is_gamma_calc:
                            if abs(one_peak.gamma_ref - one_peak.ccr_rate)>=3:
                                one_exp.deepResidueAnalisys(writer,peak_number,
                                                        error_residue_dict,
                                                        error_experiment_dict,
                                                        error_aa_dict)

            writer_row = csv.writer(csv_file,delimiter=",")
            writer_row.writerow("")
            
            self.deepExpAnalisys(writer_row,
                    error_residue_dict,
                    error_experiment_dict,
                    error_aa_dict,
                    invisible_cross_peaks,
                    aa_dict_nun)
            
            
        # print ("All info from peaks for {} rates are in file: {}".format(ccr_name, new_list), file=RaportBox)
        RaportBox.write("Information about very small CCR rates are in file: {}\n".format(new_list))
        return







class CCRClass:
    __slots__ = ['_CCR_name', '_ref_name','_trans_name','_noise','_ns', '_n_dim','_CCR_pos',
                 '_tc_vol','_rate_mult','_is_peak_uncertainty','_peaks','_is_peaklist','_other','_is_reference','r2']
    
    def __init__(self,exp_dict):
        self._CCR_name = exp_dict["type_of_CCR"]       
        self._noise = []      
        self._ns = exp_dict["NS"]         
        self._ref_name = []       # type: list[str]
        self._trans_name = []       # type: list[str]
            
        self._n_dim = 0           
        self._CCR_pos = 100
        self._tc_vol = float(exp_dict["TC"])   
        self._rate_mult = False

        #self._is_peak_uncertainty = False
        self._peaks = []               # type: list[CResidue] 
        self._is_peaklist = False
        self._other = ""
        self._is_reference = False
        self.r2 = 'nan'


        if "ref_name"in exp_dict:
            if isinstance(exp_dict["ref_name"],list):
                self._ref_name = exp_dict["ref_name"]
            else:
                self._ref_name = [exp_dict["ref_name"]]

        elif "dir_ref" in exp_dict:
            if isinstance(exp_dict["dir_ref"],list):
                self._ref_name = f'{exp_dict["dir_ref"]}_{self._CCR_name}_a'
            else:
                self._ref_name = [f'{exp_dict["dir_ref"]}_{self._CCR_name}_a']

        if "trans_name"in exp_dict:
            if isinstance(exp_dict["trans_name"],list):
                self._trans_name = exp_dict["trans_name"]
            else:
                self._trans_name = [exp_dict["trans_name"]]

        elif "dir_trans" in exp_dict:
            if isinstance(exp_dict["dir_trans"],list):
                self._trans_name = f'{exp_dict["dir_trans"]}_{self._CCR_name}_x'
            else:
                self._trans_name = [f'{exp_dict["dir_trans"]}_{self._CCR_name}_x']

        if "CCR_pos" in exp_dict:
            self._CCR_pos = deepcopy(int(exp_dict["CCR_pos"]))
        else:
            self._CCR_pos = deepcopy(int(CCR_dict[self._CCR_name]["CCR_pos"]))

        if "rate_mult" in exp_dict:
            self._rate_mult = exp_dict["rate_mult"]
        else:
            if self._CCR_name in CCR_dict:
                if "rate_mult" in CCR_dict[self._CCR_name]:
                    self._rate_mult = CCR_dict[self._CCR_name]["rate_mult"]

        if "dimension" in exp_dict:
            self._n_dim = int(exp_dict["dimension"])  
        else:
            if self._CCR_name in CCR_dict:
                if "dim" in CCR_dict[self._CCR_name]:
                    self._rate_mult = CCR_dict[self._CCR_name]["dim"]
            else:
                print(f"Error: {self._CCR_name} - dimensionality of spectrum is unknow. Please add this parameter to experiments setup file")
                sys.exit()

        if "other" in exp_dict:
            self._other = exp_dict["other"]
            
        if "noise" in exp_dict:
            self._noise = exp_dict["noise"]

        print_raport(f"""
CCR name: {self._CCR_name}
reference name: {self._ref_name}
transfer name: {self._trans_name}
n dim: {self._n_dim}
CCR pos: {self._CCR_pos}
ns: {self._ns}
Tc: {self._tc_vol}
other: {self._other}
""")
        

    def peak_list(self)->list[CResidue]:
        return self._peaks
    
    def is_peaklist(self)->bool:
        return self._is_peaklist
    
    def is_reference(self)->bool:
        return self._is_reference
    
    def CCRname(self)->str:
        return self._CCR_name
    
    def other_info(self)->str:
        return self._other
    
    def s_dim(self)->int:
        return self._n_dim
    
    def give_info_about(self,peak_number,**kwargs):  #do usuniecia
        info_dict = {}
        one_peak = self._peaks(peak_number)
        option_dict = {
            "CCR name":self._CCR_name,
            "noise":self._noise,
            "ns reference":self._ns[0],
            "ns transfer":self._ns[1],
            "dim":self._n_dim,
            "CCR pos":self._CCR_pos,
            "Tc":self._tc_vol,
            "other":self._other,
            "nucl name":self._nucl_name,
            "nucl pos":self._nucl_pos,
            "angle num":self._n_angle,
            "angle":self._angle,
            "angle pos":self._angle_pos,
            "peak pos":one_peak.peak_pos,
            "peak pos point":one_peak.peak_points_pos,
            "peak aa name":one_peak.aa_name,
            "peak aa num":one_peak.aa_number,
            "peak intens":one_peak.peak_intens,
            "peak descript":one_peak.descript,
            "peak gamma ref":one_peak.gamma_ref,
            "peak Ix_theor":one_peak.Ix_theor,
            "peak Ix_theor_Ia_ratio":one_peak.Ix_theor_Ia_ratio,
            "peak Ix_theor_Ia_ratio_without_NS":one_peak.Ix_theor_Ia_ratio_without_NS,
        }
        for one_kwarg in kwargs:
            if one_kwarg in option_dict:
                info_dict[one_kwarg] = option_dict[one_kwarg]
            else:
                print_raport(f"There is no option {one_kwarg}")
        return info_dict
    
    def Prepare_peaklist(self, seq_dict:dict,):
        for aminoacids_number in seq_dict:
            res = CResidue(aa_num = aminoacids_number,
                           aa_name = seq_dict[aminoacids_number])
            self._peaks.append(deepcopy(res))

    @staticmethod
    def read_aminoacid_number(peak_description:str) -> int:
        try:
            aminoacids_number = int(peak_description)
        except:
            aminoacids = peak_description.split("-")
            first_dim_sign_list = list(aminoacids[0])
            last_dim_sign_list = list(aminoacids[-1])
            # RaportBox.write(f"{first_dim_sign_list}, {last_dim_sign_list}")
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
            
            # RaportBox.write(f"\ndigit_list_first: {digit_list_first}, digit_list_last: {digit_list_last}\n")
            # RaportBox.write(f"aminoacids_number_first: {aminoacids_number_first}, aminoacids_number_last: {aminoacids_number_last}\n")

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
                print_raport(f"Problem with number reading in: {first_dim_sign_list} and {last_dim_sign_list}")
                aminoacids_number = "0"
            # print(f"aminoacids_number: {aminoacids_number}")
        return aminoacids_number

    def Read_peak_uncertainty(self, peaklist_directory:str, list_verson:str, version:int):
        peaklistfile = f"{peaklist_directory}/peak_lists/{list_verson}_peaks_noise.list"

        if os.path.exists(peaklistfile):
            RaportBox.write(f"\n\nLIST:{peaklistfile}\n")
            with open(peaklistfile, 'r') as peak_noise_file:  
                p_lines_a = peak_noise_file.readlines()
                for indexl, line in enumerate(p_lines_a):
                    if indexl > 0 and len(line) > 1:
                        item_a = line.split()
                        description = item_a[0]
                        aminoacids_number = self.read_aminoacid_number(description)
                        RaportBox.write(f"aminoacids_number: {aminoacids_number,description}\n")
                        aminoacids_number = int(aminoacids_number)
                        for residue in self._peaks:
                            if aminoacids_number == residue.aa_number:
                                residue.peak_uncertainty[version]=(deepcopy(float(item_a[-1])))
                                RaportBox.write(f"peak_uncertainty = {float(item_a[-1])}\n")
            # self._is_peak_uncertainty = True
        else:
            # self._is_peak_uncertainty = False
            for residue in self._peaks:
                residue.peak_uncertainty=[self._noise[0],self._noise[1]]
            print_raport (f"There is no file: {peaklistfile}!, so average noise value will be used\n")
    
    @abstractmethod
    def Read_peaklist(self, peaklist_directory, seq_dict:dict, points_mode=False):
        pass

    @abstractmethod
    def read_peak_files(self,peaklist_directory:str,seq_dict:dict):
        pass
        
    @abstractmethod
    def calc_uncertainty_value(self,peak_number):
        pass
        
    @abstractmethod
    def calc_ccr_rate(self):
        pass
    
    @abstractmethod
    def calc_theor_Ix(self):
        pass


    @staticmethod
    def calc_distance(peak1:list[int],peak2:list[int])->float:
        sum_of_squares = 0
        for i in range(len(peak1)):
            # print (peak1, "and", peak2)
            sum_of_squares += (peak1[i]-peak2[i])**2
        dis =  sum_of_squares**(1/len(peak1))
        return dis
         
    def check_overlap(self, peak1:CResidue, peak2:CResidue, version=0):
        peaks_distance = self.calc_distance(peak1.peak_pos_points[version], peak2.peak_pos_points[version])
        if peaks_distance <= 4:
            RaportBox.write(f"\nCheck overlap for {peak1.descript}: {peak1.peak_pos_points[version]} and {peak2.descript}:{peak2.peak_pos_points[version]} - distance: {peaks_distance:.3f}")
            if peaks_distance < 2:
                peak1.is_overlap = True
                peak2.is_overlap = True
                peak1.overlap_peaks[f"{peak2.aa_number}_{version}"]=peaks_distance
                peak2.overlap_peaks[f"{peak1.aa_number}_{version}"]=peaks_distance
                RaportBox.write(f"\nPeaks: {peak1.descript} and {peak2.descript} marked as overlapping")
            


    def check_overlap_for_all_peaks(self):
        for indexp1, peak1 in enumerate(self._peaks):
            if peak1.is_peak:
                for indexp2, peak2 in enumerate(self._peaks):
                    if indexp2 < indexp1:
                        if peak2.is_peak:
                            for i in range(len(peak1.peak_pos)):
                                self.check_overlap(peak1,peak2,version=i)

    def Add_ref_gamma(self,ref_seq_num,ref_CCR):
        # print ("ref_seq_num",ref_seq_num)
        for one_peak in self._peaks:
            for indexr, r_seq_num in enumerate(ref_seq_num):
                # print ("aa_number: {}, seq_number: {}".format(one_peak.aa_number,r_seq_num))
                if int(one_peak.aa_number) == int(r_seq_num) and ref_CCR[indexr] != 'nan':
                    one_peak.is_gamma_calc = True
                    one_peak.gamma_ref = float(ref_CCR[indexr])
                    self._is_reference = True
                    # print ("gamma_ref",one_peak.gamma_ref)
    
    def check_if_fatal_error(self):
        min_max_value = [+100.0,-100.0] 
        for one_peak in self._peaks:
            if one_peak.is_gamma_calc and one_peak.is_ccr_rate:
                min_max_value = check_if_min_max(one_peak.gamma_ref,one_peak.ccr_rate,min_max_value)
        
        peak_range = math.sqrt((min_max_value[1] - min_max_value[0])**2)
        for one_peak in self._peaks:    
            if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.ccrrate_calculation_error == False:
                dist = math.sqrt((one_peak.ccr_rate - one_peak.gamma_ref )**2)
                if dist > peak_range*0.2:
                    one_peak.fatal_error = True

    def Additional_text(self)->str:
        if self._other=="":
            add_text=""
        else: add_text="_"+self._other
        return add_text

    def plot_gamma_gamma(self, transparent_plot=False, add="",add2="",style='none'):
        if add2 == "":
            add2 = self._other
        min_max_value = [+100.0,-100.0]         #minimal and maximal value of CCR rate or reference gamma - it is necessary to make plot  
        gamma_calculated = []
        gamma_experimental = []
        gamma_calc_error = []            # ERROR GAMMA INFO, 9.11.2023
        gamma_differences = []
        for one_peak in self._peaks:
            if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
                # print(one_peak.gamma_ref,one_peak.ccr_rate)
                gamma_calculated.append(deepcopy(one_peak.gamma_ref))
                gamma_experimental.append(deepcopy(one_peak.ccr_rate))
                gamma_calc_error.append(deepcopy(one_peak.ccrrate_error_value))            # ERROR GAMMA INFO, 9.11.2023
                gamma_differences.append(deepcopy(one_peak.ccr_rate-one_peak.gamma_ref))
                min_max_value = check_if_min_max(one_peak.gamma_ref,one_peak.ccr_rate,min_max_value)

        plt.rcParams['font.size'] = '14'

        wilcoxon_rank = stats.wilcoxon(gamma_calculated,gamma_experimental)[1] # results: statistic, pvalue, zstatistic  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
        tStudent_2 = stats.ttest_ind(gamma_calculated,gamma_experimental)[0]    # results: statistic, pvalue, df  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
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
       
        # print_raport(f"\n {self._CCR_name}")
        RaportBox.write(label)
        if style == "publication":
            plt.title(f'{CCRname2PrettyRateNamePLT(self.CCRname())}', fontsize=10)
        elif style == "presentation":
            plt.title(f'{CCRname2PrettyRateNamePLT(self.CCRname())}', fontsize=10)
        else:
            plt.title(self._CCR_name+add)
            # print(f"gamma_calculated:{gamma_calculated} \ngamma_experimental: {gamma_experimental}\ngamma_calc_error: {gamma_calc_error}")
            linear_expretion,linear_r2, linear_slope, linear_intercept, matching_factor_for_y_x_LR= LRegression_expresion(gamma_calculated, gamma_experimental)
            weighted_reg_dict = WeightedLRegression_expresion_by_hand(gamma_calculated, gamma_experimental,gamma_calc_error)
            self.r2 = weighted_reg_dict["r2"]
            plt.axline([0,linear_intercept],slope=linear_slope, linestyle=(0, (5, 5)), 
                       linewidth=1.5, color='red', 
                       label=f'l: {linear_expretion}, r2 = {linear_r2:.3f}, \nmatching factor for y=x: {matching_factor_for_y_x_LR:.1f}')
            plt.axline([0,weighted_reg_dict["intercept"]],slope=weighted_reg_dict["slope"], 
                       linestyle=(0, (5, 5)), linewidth=1.5, color='purple', 
                       label=f'lw2: {weighted_reg_dict["equation"]}, \nmatching factor for y=x: {weighted_reg_dict["factor_a"]:.1f}')
            plt.legend(fontsize="small",bbox_to_anchor=(1.05, 1),
                            loc='upper left', borderaxespad=0.)

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
        
        # To load the display window
        # plt.savefig("{}/wykresy2/{}_exp_vs_calc.png".format(Dir, tfn[0].name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
        plt.savefig("{}/{}_exp_vs_calc{}.png".format(file_directory, self._CCR_name,add2), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
        # plt.show()
        plt.clf()
        plt.close()
    



    def plot_gamma_gamma_with_fatal(self, transparent_plot=False, add="",add2="",style='none'):
        if add2 == "":
            add2 = self._other
        min_max_value = [+100.0,-100.0]         #minimal and maximal value of CCR rate or reference gamma - it is necessary to make plot  
        gamma_calculated = []
        gamma_experimental = []
        gamma_calc_error = []            # ERROR GAMMA INFO, 9.11.2023
        gamma_differences = []
        fatal_error_point_calc = []
        fatal_error_point_exp = []
        fatal_error_point_err = []
        for one_peak in self._peaks:
            if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G":
                # print(one_peak.gamma_ref,one_peak.ccr_rate)
                if one_peak.fatal_error:
                    fatal_error_point_calc.append(deepcopy(one_peak.gamma_ref))
                    fatal_error_point_exp.append(deepcopy(one_peak.ccr_rate))
                    fatal_error_point_err.append(deepcopy(one_peak.ccrrate_error_value))
                else:
                    gamma_calculated.append(deepcopy(one_peak.gamma_ref))
                    gamma_experimental.append(deepcopy(one_peak.ccr_rate))
                    gamma_calc_error.append(deepcopy(one_peak.ccrrate_error_value))            # ERROR GAMMA INFO, 9.11.2023
                    gamma_differences.append(deepcopy(one_peak.ccr_rate-one_peak.gamma_ref))
                min_max_value = check_if_min_max(one_peak.gamma_ref,one_peak.ccr_rate,min_max_value)

        plt.rcParams['font.size'] = '14'

        wilcoxon_rank = stats.wilcoxon(gamma_calculated,gamma_experimental)[1] # results: statistic, pvalue, zstatistic  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
        tStudent_2 = stats.ttest_ind(gamma_calculated,gamma_experimental)[0]    # results: statistic, pvalue, df  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
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
        #fatal errors points:
        if len(fatal_error_point_calc)>0 and len(fatal_error_point_calc)>0 and len(fatal_error_point_calc)>0:
            plt.scatter(fatal_error_point_calc,fatal_error_point_exp,s=3, color='#252525ff') #
            plt.errorbar(fatal_error_point_calc, fatal_error_point_exp, 
                    yerr=fatal_error_point_err, 
                    fmt='none', color='#252525ff') #

        # print_raport(f"\n {self._CCR_name}")
        RaportBox.write(label)
        if style == "publication":
            plt.title(f'{CCRname2PrettyRateNamePLT(self.CCRname())}', fontsize=10)
        elif style == "presentation":
            plt.title(f'{CCRname2PrettyRateNamePLT(self.CCRname())}', fontsize=10)
        else:
            plt.title(self._CCR_name+add)
            # print(f"gamma_calculated:{gamma_calculated} \ngamma_experimental: {gamma_experimental}\ngamma_calc_error: {gamma_calc_error}")
            linear_expretion,linear_r2, linear_slope, linear_intercept, matching_factor_for_y_x_LR= LRegression_expresion(gamma_calculated, gamma_experimental)
            weighted_reg_dict = WeightedLRegression_expresion_by_hand(gamma_calculated, gamma_experimental,gamma_calc_error)
            self.r2 = weighted_reg_dict["r2"]
            plt.axline([0,linear_intercept],slope=linear_slope, linestyle=(0, (5, 5)), 
                       linewidth=1.5, color='red', 
                       label=f'l: {linear_expretion}, r2 = {linear_r2:.3f}, \nmatching factor for y=x: {matching_factor_for_y_x_LR:.1f}')
            plt.axline([0,weighted_reg_dict["intercept"]],slope=weighted_reg_dict["slope"], 
                       linestyle=(0, (5, 5)), linewidth=1.5, color='purple', 
                       label=f'lw2: {weighted_reg_dict["equation"]}, \nmatching factor for y=x: {weighted_reg_dict["factor_a"]:.1f}')
            plt.legend(fontsize="small",bbox_to_anchor=(1.05, 1),
                            loc='upper left', borderaxespad=0.)

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
        
        # To load the display window
        # plt.savefig("{}/wykresy2/{}_exp_vs_calc.png".format(Dir, tfn[0].name), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
        plt.savefig("{}/{}_exp_vs_calc{}.png".format(file_directory, self._CCR_name,add2), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
        # plt.show()
        plt.clf()
        plt.close()
    
    def plot_error_histogram(self,transparent_plot=False,add="",add2=""): # ERROR GAMMA INFO, 9.11.2023
        gamma_experimental = []
        gamma_calc_error = []            
        peak_intens_auto = []         
        peak_intens_cross = []
        for one_peak in self._peaks:
            if one_peak.is_ccr_rate and one_peak.is_gamma_calc and one_peak.aa_name!="G" and one_peak.is_peak:
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
        axs[0].set_title("{} {}\n Experimental \u0393 vs Uncertainty value".format(self._CCR_name,add))
        axs[0].set_ylabel('Uncertainty value')
        axs[0].set_xlabel('Experimental \u0393, $s^{-1}$')

        axs[1].scatter(peak_intens_auto, gamma_calc_error, s=5, color='#0066ffff')
        axs[1].set_title("Peak intens in reference vs Uncertainty value")
        axs[1].set_ylabel('Uncertainty value')
        axs[1].set_xlabel('Peak intens in reference')

        axs[2].scatter(peak_intens_cross, gamma_calc_error, s=5, color='#0066ffff')
        axs[2].set_title("Peak intens in transfer vs Uncertainty value")
        axs[2].set_ylabel('Uncertainty value')
        axs[2].set_xlabel('Peak intens in transfer')

        # plt.title(ccr_name+add)
        plt.savefig("{}Hist_{}_exp_vs_Uncertainty{}.png".format(RaportDir, self._CCR_name,add2), bbox_inches="tight", pad_inches=0.3, transparent=transparent_plot)
        # plt.show()
        plt.clf()
        plt.close()

    def deepResidueAnalisys(self, writer, peak_number:int,
                    error_residue_dict:dict,
                    error_experiment_dict:dict,
                    error_aa_dict:dict):
        """ For particular peak we write all information about it and append basics info inf dicts
        
        """
        CCR_name =  self._CCR_name+self.Additional_text()
        one_peak = self._peaks[peak_number]
        peak_descr = str(one_peak.aa_number)+one_peak.aa_name
        one_row = {}
        error_experiment_dict[CCR_name]+=1
        if one_peak.aa_name in error_aa_dict:
            error_aa_dict[one_peak.aa_name]+=1
        else:
            error_aa_dict[one_peak.aa_name]=1
        one_row["CCR name"] = CCR_name
        one_row["AA"] = peak_descr

        if one_peak.is_peak == True and one_peak.is_ccr_rate == True:
            if self._CCR_pos == 0:
                for i in range(self._n_dim):
                    one_row['w{}'.format(i+1)] = '{:.3f}'.format(one_peak.peak_pos[0][i])
                one_row["intensity in reference"] = one_peak.peak_intens[0]
                one_row["intensity in transfer"] = one_peak.peak_intens[1]
                one_row["Noise in reference"] = one_peak.peak_uncertainty[0]
                one_row["Noise in transfer"] = one_peak.peak_uncertainty[1]
                                
                if one_peak.is_peak and one_peak.is_theor_values:
                    one_row['Ix Theor'] = "{:.0f}".format(one_peak.Ix_theor)
                    one_row['Ix/Ia'] = "{:.4f}".format(one_peak.Ix_theor_Ia_ratio)
                    one_row['Ix/Ia without NS'] = "{:.4f}".format(one_peak.Ix_theor_Ia_ratio_without_NS)
                    
            if one_peak.ccrrate_calculation_error == False:
                one_row["CCR rate"] = '{:.4f}'.format(one_peak.ccr_rate)
            elif one_peak.to_check:
                one_row["Comments"] = "Check peak position!"
            else:
                one_row["Comments"] = "peak intens error"

        elif one_peak.is_peak == False and one_peak.is_ccr_rate == True: 
            if one_peak.aa_name == "G":
                one_row["Comments"] = "glycine - no CCR rate"
            elif one_peak.to_check:
                one_row["Comments"] = "Check peak position!"
            else:
                one_row["CCR rate"] = '{:.4f}'.format(one_peak.ccr_rate)
                one_row["Comments"] = "Info from other peak"
        elif one_peak.is_overlap == True:
            one_row["Comments"] = "This peak overlap with: {}".format(one_peak.overlap_peaks)
        elif one_peak.is_peak == True and one_peak.is_ccr_rate == False:
            for i in range(self._n_dim):
                one_row['w{}'.format(i+1)] = '{:.3f}'.format(one_peak.peak_pos[0][i])
            if one_peak.to_check:
                        one_row["Comments"] = "Check peak position!"
            # if one_peak.is_theor_values:
            #     one_row['Ix Theor'] = "{:.0f}".format(one_peak.Ix_theor)
            #     one_row['Ix/Ia'] = "{:.4f}".format(one_peak.Ix_theor_Ia_ratio)
            #     one_row['Ix/Ia without NS'] = "{:.4f}".format(one_peak.Ix_theor_Ia_ratio_without_NS)
        else:
            one_row["Comments"] = "No visible peak"
        
        if one_peak.ccrrate_error_value != 0.0:
            one_row['Uncertainty'] = '{:.4f}'.format(one_peak.ccrrate_error_value)            # ERROR GAMMA INFO, 9.11.2023
        

        if one_peak.is_gamma_calc and one_peak.is_ccr_rate:
            one_row['Reference Gamma'] = "{:.4f}".format(one_peak.gamma_ref)

        
        
        writer.writerow(one_row)

        if self._CCR_pos == -1:
            next_peak = self._peaks[peak_number+1]
            one_row = {}
            peak_descr = str(next_peak.aa_number)+next_peak.aa_name
            one_row["CCR name"] = CCR_name
            one_row["AA"] = peak_descr
            for i in range(self._n_dim):
                one_row['w{}'.format(i+1)] = '{:.3f}'.format(next_peak.peak_pos[0][i])
            one_row["intensity in reference"] = next_peak.peak_intens[0]
            one_row["intensity in transfer"] = next_peak.peak_intens[1]
            one_row["Noise in reference"] = self._noise[0]
            one_row["Noise in transfer"] = self._noise[1]

            if next_peak.is_peak and next_peak.is_theor_values:
                one_row['Ix Theor'] = "{:.0f}".format(next_peak.Ix_theor)
                one_row['Ix/Ia'] = "{:.4f}".format(next_peak.Ix_theor_Ia_ratio)
                one_row['Ix/Ia without NS'] = "{:.4f}".format(next_peak.Ix_theor_Ia_ratio_without_NS)
            
            writer.writerow(one_row)

            if peak_descr in error_residue_dict:
                error_residue_dict[peak_descr].append(deepcopy(CCR_name+" (-1 position)"))
            else:
                error_residue_dict[peak_descr] = [CCR_name+" (-1 position)"]

        elif self._CCR_pos == 0:
            if peak_descr in error_residue_dict:
                error_residue_dict[peak_descr].append(deepcopy(CCR_name))
            else:
                error_residue_dict[peak_descr] = [CCR_name]

            

    def deepVisibilityAnalisis(self,one_peak:CResidue,invisible_cross_peaks):
        """Checking if peak if visible in transfer version of experiment
        """
        CCR_name =  self._CCR_name+self.Additional_text()
        noise_level = self._noise
        peak_descr = str(one_peak.aa_number)+one_peak.aa_name
        if one_peak.is_peak:
            if abs(one_peak.peak_intens[1])<noise_level[1]*10:
                invisible_cross_peaks[CCR_name][0]+=1
                RaportBox.write("{} - {}, noise: {} intens: {} \n".format(invisible_cross_peaks[CCR_name][0],peak_descr,noise_level[1],one_peak.peak_intens[1]))
            if abs(one_peak.peak_intens[1])<noise_level[1]:
                invisible_cross_peaks[CCR_name][1]+=1
    
    def WriteCCRRate_small(self):
        add = self.Additional_text()
        new_list = "{}{}{}.csv".format(RaportDir,CCRname2Ratename(self._CCR_name),add)
        with open(new_list, mode='w', newline='') as csv_file:
            headers = ['AA','CCR rate']
            writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
            for one_peak in self._peaks:
                one_row = {}
                one_row["AA"] = one_peak.aa_number
                if one_peak.is_ccr_rate and not one_peak.is_overlap and not one_peak.ccrrate_calculation_error and one_peak.aa_name != "G":
                    one_row["CCR rate"] = '{:.4f}'.format(one_peak.ccr_rate)
                else:
                    one_row["CCR rate"] = "nan"
                writer.writerow(one_row)
        # print ("CCR rate for CCR effect between {} rates are in file: {}".format(CCRname2Ratename(ccr_name), new_list), file=RaportBox)
        RaportBox.write("CCR rate for CCR effect between {} rates are in file: {}\n".format(CCRname2Ratename(self._CCR_name), new_list))
        return

    def WriteCCRRate_all_info(self):
        add = self.Additional_text()
        new_list = "{}{}{}_CCRrates.list".format(RaportDir,self._CCR_name,add)
        # print ("new_peak_list", new_list, file=RaportBox)
        max_lenth_discrip = 0
        max_lenth_intens = 0

        for p in self._peaks:
            if p.is_peak:
                pd = str(p.aa_number)+p.aa_name
                if len(pd)>max_lenth_discrip:
                    max_lenth_discrip=len(pd)
                if len(str(p.peak_intens[0]))>max_lenth_intens:
                    max_lenth_intens=len(str(p.peak_intens[0]))
        
        with open(new_list, 'w') as listfile:
            sentence_len_7_dim = int(7*self._n_dim)
            print ("\tauto list name = {}\tcross list name = {}\tTc = {}".format(self._ref_name[0], self._trans_name[0], self._tc_vol), file=listfile) 
            print ("{0}\t{1:4}\t{2:5}\t{3:5}\t CCR rate \t CCR Uncertainty \t CCR theor \t Ix theor \t Ix/Ia \t Comments".format("AA",
                                                                                                                    "peak position",
                                                                                                                    "intensity in reference",
                                                                                                                    "intensity in transfer",
                                                                                                                    sentence_len_7_dim,
                                                                                                                    max_lenth_intens)
                                                                                                                    , file=listfile)
            for one_peak in self._peaks:
                dict_line={}
                peak_descr = str(one_peak.aa_number)+one_peak.aa_name
                dict_line['AA'] = "\n{:{sentence_len}}".format(peak_descr, sentence_len=max_lenth_discrip)
                # print ("{:{sentence_len}}".format(peak_descr, sentence_len=max_lenth_discrip,), end="\t", file=listfile)
                if one_peak.is_peak == True:
                    peak_pos = ""
                    for i in range(self._n_dim):
                        peak_pos += "{:.3f}\t".format(one_peak.peak_pos[0][i])
                    dict_line['peak_pos'] = peak_pos
                    dict_line['Ia'] = "{0:{sentence_len}}".format("{:.2e}".format(one_peak.peak_intens[0]), sentence_len=max_lenth_intens) 
                    dict_line['Ix'] = "{0:{sentence_len}}".format("{:.2e}".format(one_peak.peak_intens[1]), sentence_len=max_lenth_intens) 
                else:
                    if one_peak.is_ccr_rate == True:
                        dict_line['peak_pos'] = "{0:{sentence_len}}\t\t".format("Info from other peak", sentence_len=7*self._n_dim+max_lenth_intens*2)
                    else:
                        dict_line['peak_pos'] = "{0:{sentence_len}}\t\t".format("No visible peak", sentence_len=7*self._n_dim+max_lenth_intens*2) 
                

                if one_peak.aa_name == "G":
                    dict_line['ccr_rate'] = "{:{sentence_len}}".format("glycine", sentence_len=18)   #21 digit, glycine - no CCR rate
                # elif one_peak.is_overlap == True:
                #     dict_line['ccr_rate'] = "This peak overlap with: {}".format(one_peak.overlap_peaks)
                elif one_peak.ccrrate_calculation_error == True:
                    dict_line['ccr_rate'] = "peak intens error: Ix/Ia = {:8.4f}".format(one_peak.is_ccr_rate, sentence_len=18) 
                # elif one_peak.is_peak == True and one_peak.is_ccr_rate == False:
                #     dict_line['ccr_rate'] = "{:{sentence_len}}".format("Info to other peak", sentence_len=18)
                elif one_peak.ccrrate_calculation_error == False and one_peak.is_ccr_rate and one_peak.is_overlap == False:
                    dict_line['ccr_rate'] = '{:18.4f}'.format(one_peak.ccr_rate, sentence_len=18)
                else:
                    dict_line['ccr_rate'] = "         -          "
                

                if one_peak.ccrrate_error_value != 0.0:
                    dict_line['Uncertainty'] = '{:8.4f}'.format(one_peak.ccrrate_error_value, sentence_len=10)            # ERROR GAMMA INFO, 9.11.2023
                else:
                    dict_line['Uncertainty'] = " "*10  
                

                if one_peak.is_gamma_calc and one_peak.is_ccr_rate:
                    dict_line['ccr_rate_theor'] = "{:8.4f}".format(one_peak.gamma_ref, sentence_len=10)
                else:
                    dict_line['ccr_rate_theor'] =  dict_line['Ix_theor'] = dict_line['Ix_Ia'] = " "*10  

                if one_peak.is_theor_values:
                    dict_line['Ix_theor'] = "{:8.2e}".format(one_peak.Ix_theor, sentence_len=10)
                    dict_line['Ix_Ia'] = "{:8.2f}".format(one_peak.Ix_theor_Ia_ratio, sentence_len=10)
                else:
                    dict_line['Ix_theor'] =  dict_line['Ix_theor'] = dict_line['Ix_Ia'] = " "*10  
                    dict_line['Ix_Ia'] =  dict_line['Ix_theor'] = dict_line['Ix_Ia'] = " "*10  

                dict_line["Comments"] = ""

                if one_peak.to_check:
                    dict_line["Comments"] += "Check peak position! "
                if one_peak.is_overlap == True:
                    dict_line['ccr_rate'] += "This peak overlap with: {} ".format(one_peak.overlap_peaks)
                if one_peak.is_peak == False and one_peak.is_ccr_rate:
                    dict_line['ccr_rate'] = "{:{sentence_len}}".format("CCR rate value from peak: {} ".format(), 
                                                                       sentence_len=18)


                for item in dict_line:
                    print (dict_line[item], end="\t", file=listfile)
                
                
        # print ("All info from peaks for {} rates are in file: {}".format(ccr_name, new_list), file=RaportBox)
        RaportBox.write("All info from peaks for {} rates are in file: {}\n".format(self._CCR_name, new_list))
        return

    def WriteCCRRate_all_info_CSV(self):
        add = self.Additional_text()
        new_list = "{}{}{}_CCRrates.csv".format(RaportDir,self._CCR_name,add)
        max_lenth_discrip = 0
        max_lenth_intens = 0
        for p in self._peaks:
            if p.is_peak:
                pd = str(p.aa_number)+p.aa_name
                if len(pd)>max_lenth_discrip:
                    max_lenth_discrip=len(pd)
                if len(str(p.peak_intens[0]))>max_lenth_intens:
                    max_lenth_intens=len(str(p.peak_intens[0]))
        with open(new_list, mode='w', newline='') as csv_file:
            # headers = ['AA','peak position','intensity in reference','intensity in transfer','CCR rate','Comments']
            headers = ['AA']
            for i in range(1,self._n_dim+1):
                headers.append(deepcopy('w{}'.format(i)))
            headers.extend(deepcopy(['intensity in reference','intensity in transfer','CCR rate','Uncertainty','Comments','Reference Gamma', 'Ix Theor', 'Ix/Ia', 'Ix/Ia without NS']))
            writer_row = csv.writer(csv_file,delimiter=",")
            writer_row.writerow(['reference list name =',self._ref_name[0]])
            writer_row.writerow(['transfer list name =',self._trans_name[0]])
            writer_row.writerow(['Tc =', self._tc_vol])
            writer_row.writerow(['---','---','---','---','---','---','---','---','---'])
            # writer_row.writerow(['AA','{:^{sentence_len1}}'.format('peak position',sentence_len1=7*s_dim),'intensity in reference','intensity in transfer','CCR rate','Comments'])
            writer = csv.DictWriter(csv_file, fieldnames=headers, delimiter=",")
            writer.writeheader()
            for one_peak in self._peaks:
                one_row = {}
                peak_descr = str(one_peak.aa_number)+one_peak.aa_name
                one_row["AA"] = peak_descr
                if one_peak.is_peak == True and one_peak.is_ccr_rate == True:
                    for i in range(self._n_dim):
                        one_row['w{}'.format(i+1)] = '{:.3f}'.format(one_peak.peak_pos[0][i])
                    if one_peak.aa_name == "G":
                        one_row["Comments"] = "glycine - no CCR rate"
                    else:
                        one_row["intensity in reference"] = one_peak.peak_intens[0]
                        one_row["intensity in transfer"] = one_peak.peak_intens[1]
                        if one_peak.ccrrate_calculation_error == False:
                            one_row["CCR rate"] = '{:.4f}'.format(one_peak.ccr_rate)
                        elif one_peak.to_check:
                            one_row["Comments"] = "Check peak position!"
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
                    for i in range(self._n_dim):
                        one_row['w{}'.format(i+1)] = '{:.3f}'.format(one_peak.peak_pos[0][i])
                    if one_peak.to_check:
                        one_row["Comments"] = "Check peak position!"
                else:
                    one_row["Comments"] = "No visible peak"
                
                if one_peak.ccrrate_error_value != 0.0:
                    one_row['Uncertainty'] = '{:.4f}'.format(one_peak.ccrrate_error_value)            # ERROR GAMMA INFO, 9.11.2023
                

                if one_peak.is_gamma_calc and one_peak.is_ccr_rate:
                    one_row['Reference Gamma'] = "{:.4f}".format(one_peak.gamma_ref)
                if one_peak.is_theor_values:
                    one_row['Ix Theor'] = "{:.2e}".format(one_peak.Ix_theor)
                    one_row['Ix/Ia'] = "{:.2f}".format(one_peak.Ix_theor_Ia_ratio)
                    one_row['Ix/Ia without NS'] = "{:.2f}".format(one_peak.Ix_theor_Ia_ratio_without_NS)
                writer.writerow(one_row)
        # print ("All info from peaks for {} rates are in file: {}".format(ccr_name, new_list), file=RaportBox)
        RaportBox.write("All info from peaks for {} rates are in file: {}\n".format(self._CCR_name, new_list))
        return
    
    @abstractmethod
    def WriteCCRRate_all_info_new(self, file_type:str):
        """ Methods for writing info for one CCR rate.

        file_type: "csv" or "txt"
        """
        pass

class CCR_normal(CCRClass):
    
    def __init__(self,exp_dict):
        super().__init__(exp_dict)

    def Read_peaklist(self, peaklist_directory, points_mode=False):
        NameFlag = [False, False]
        if points_mode == False:
            list_of_names_ends = [".list","_new_ppm.list"]
        else:
            list_of_names_ends = ["_points.list","_new_points.list"]
        peak_list_basic_names = [self._ref_name[0], self._trans_name[0]]
        peak_list_names = ["None","None"]
        for indexlv, list_verson in enumerate(peak_list_basic_names):
            for i, l in enumerate(list_of_names_ends):
                peaklistfile = f"{peaklist_directory}/{list_verson}{l}"
                try:
                    f=open(peaklistfile)
                    NameFlag[indexlv] = True
                    peak_list_names[indexlv] = peaklistfile
                    f.close
                except FileNotFoundError:
                    RaportBox.write("There is no such file or directory: {}\n".format(peaklistfile))
        
        if NameFlag[0] == True and NameFlag[1] == True:
            RaportBox.write("\n\nLISTA:{}\n".format(peak_list_names))
            with open(peak_list_names[0], 'r') as pl_a, open(peak_list_names[1], 'r') as pl_x:  
                self._is_peaklist = True
                p_lines_a = pl_a.readlines()
                p_lines_x = pl_x.readlines()
                if points_mode:
                # if points_mode and "Average noise level" in p_lines_a[1] and p_lines_x[1]:
                    self._noise = [float(p_lines_a[1].split("=")[1]), float(p_lines_x[1].split("=")[1])]
                    print_raport(f"{self._CCR_name}: auto_noise = {self._noise[0]:.2e}, cross_noise = {self._noise[1]:.2e}")
                    # RaportBox.write(f"auto_noise = {self._noise[0]}, cross_noise = {self._noise[1]}")
                for indexl, line in enumerate(p_lines_a):
                    if indexl > 1 and len(line) > 1:
                        item_a = line.split()
                        item_x = p_lines_x[indexl].split()
                        # Reading description
                        if item_a[0] != item_x[0]:
                            print ("Error: reference ({}) and transfer ({}) peak lists are not compatible:\n{} ? {}".format(peak_list_names[0], peak_list_names[1], item_a[0], item_x[0]))
                            sys.exit()
                        description = item_a[0]
                        aminoacids_number = self.read_aminoacid_number(description)
                        RaportBox.write(f"aminoacids_number: {aminoacids_number,description}\n")
                        aminoacids_number = int(aminoacids_number)
                        for res in self._peaks: 
                            if aminoacids_number == res.aa_number:
                                res.is_peak = True
                                res.descript = deepcopy(description)
                                if points_mode == True:  # for points mode
                                    for i in range(1,self._n_dim+1):
                                        res.peak_pos_points[0].append(deepcopy(int(item_a[i])))
                                        res.peak_pos_points[1].append(deepcopy(int(item_x[i])))
                                        # print ("residues append --->", residue.peak_pos_points)
                                else: 
                                    for i in range(1,self._n_dim+1):
                                        res.peak_pos[0].append(deepcopy(float(item_a[i])))
                                        res.peak_pos[1].append(deepcopy(float(item_x[i])))
                        # Reading peak intensity
                                if len(item_a)>self._n_dim+1:
                                    RaportBox.write("reference = {}, transfer = {}\n".format(item_a[self._n_dim+1],item_x[self._n_dim+1]))
                                    res.peak_intens = [ float(item_a[self._n_dim+1]) , float(item_x[self._n_dim+1]) ]

                                if len(item_a)>self._n_dim+2:
                                    if item_a[self._n_dim+2]=="check_pos" or item_a[self._n_dim+2]=="to_check":
                                        res.to_check = True
                                if len(item_x)>self._n_dim+2:
                                    if item_x[self._n_dim+2]=="check_pos" or item_x[self._n_dim+2]=="to_check":
                                        res.to_check = True
                                        
        else:
            print_raport ("Missing file for {} and {}! {}\n".format(self._ref_name[0], self._trans_name[0], peak_list_names))
            self._is_peaklist = False
    
    def read_peak_files(self,peaklist_directory:str,seq_dict:dict):
        self.Prepare_peaklist(seq_dict)
        for res in self._peaks: 
            res.peak_uncertainty = [0.0, 0.0]
            res.peak_pos = [[],[]]
            res.peak_pos_points = [[],[]]
        self.Read_peaklist(peaklist_directory, points_mode=True)
        self.Read_peaklist(peaklist_directory, points_mode=False)
        self.Read_peak_uncertainty(peaklist_directory, self._ref_name[0], version=0)
        self.Read_peak_uncertainty(peaklist_directory, self._trans_name[0], version=1)
        
    def calc_uncertainty_value(self,peak_number):
        Ia = self._peaks[peak_number].peak_intens[0]
        Ix = self._peaks[peak_number].peak_intens[1]
        NSa = self._ns[0]
        NSx = self._ns[1]
        noiseIa = self._noise[0]
        noiseIx = self._noise[1]
        uncertainty_Ia = self._peaks[peak_number].peak_uncertainty[0]
        uncertainty_Ix = self._peaks[peak_number].peak_uncertainty[1]
        RaportBox.write(f"\n{self._peaks[peak_number].descript}: Peak_uncertainty for CalcErrorValue: {uncertainty_Ia}, {uncertainty_Ix}")

        squareDerivativeArcTanh = math.pow(1/(1-math.pow((Ix*NSa)/(Ia*NSx),2)),2)
        sumOfSquaresOfPeakHightDeviation = math.pow(noiseIx/Ia,2) + math.pow(noiseIa*Ix/math.pow(Ia,2),2)
        squareError = math.pow(1/self._tc_vol,2) * squareDerivativeArcTanh * math.pow(NSa/NSx,2) * sumOfSquaresOfPeakHightDeviation
        self._peaks[peak_number+self._CCR_pos].ccrrate_error_value = math.sqrt(squareError)

        RaportBox.write("\nCalculating gamma uncertainty for peak number: {} \n squareDerivativeArcTanh = {} \n sumOfSquaresOfPeakHightDeviation = {}\n squareError = {}\n \t\t--> {}\n".format(peak_number,
                                                                                                                                                        squareDerivativeArcTanh,
                                                                                                                                                        sumOfSquaresOfPeakHightDeviation,
                                                                                                                                                        squareError,
                                                                                                                                                        math.sqrt(squareError)))
    
    def calc_ccr_rate(self):
        for indexp, peak in enumerate(self._peaks):
            if not peak.is_overlap and peak.is_peak:
                res_num = indexp+self._CCR_pos
                other_peak = self._peaks[res_num]
                other_peak.is_ccr_rate = True
                self.calc_uncertainty_value(indexp)
                try:
                    # print (peak.peak_intens[1], self.ns[0], peak.peak_intens[0],self.ns[1], "---->", (peak.peak_intens[1]*self.ns[0])/(peak.peak_intens[0]*self.ns[1]))
                    if self._rate_mult:
                        ccr_rate_vol = -atanh((peak.peak_intens[1]*self._ns[0])/(peak.peak_intens[0]*self._ns[1]))/self._tc_vol
                    else:
                        ccr_rate_vol = atanh((peak.peak_intens[1]*self._ns[0])/(peak.peak_intens[0]*self._ns[1]))/self._tc_vol
                except:
                    ccr_rate_vol = (peak.peak_intens[1]*self._ns[0])/(peak.peak_intens[0]*self._ns[1]) #"atanh(x) - x shoudl be beetween -1 and 1"
                    peak.ccrrate_calculation_error = True
                other_peak.ccr_rate = ccr_rate_vol
                
        return

    def calc_theor_Ix(self):
        for indexp, one_peak in enumerate(self._peaks):
            res_num = indexp-self._CCR_pos
            if len(self._peaks)> res_num >= 0 and one_peak.is_ccr_rate:
                other_peak = self._peaks[res_num]
                residue = str(other_peak.aa_number)+other_peak.aa_name
                if one_peak.is_gamma_calc and other_peak.is_peak:
                    # residue = str(other_peak.aa_number)+other_peak.aa_name
                    if one_peak.to_check:
                        print_raport(f"{residue} - Check peak position!!!")
                    else: print_raport(residue)

                    auto_peak_intensity = other_peak.peak_intens[0]
                    other_peak.Ix_theor = float(math.tanh(one_peak.gamma_ref*self._tc_vol)*auto_peak_intensity*self._ns[1]/self._ns[0])
                    if self._CCR_name == "CCR_5" or self._CCR_name == "CCR_6":
                        other_peak.Ix_theor = -1 * other_peak.Ix_theor
                    other_peak.Ix_theor_Ia_ratio = other_peak.Ix_theor/auto_peak_intensity
                    other_peak.Ix_theor_Ia_ratio_without_NS = other_peak.Ix_theor_Ia_ratio*self._ns[0]/self._ns[1]
                    if other_peak.Ix_theor != 'nan' and other_peak.Ix_theor_Ia_ratio != 'nan' and other_peak.Ix_theor_Ia_ratio_without_NS != 'nan':
                        other_peak.is_theor_values = True
                elif one_peak.is_gamma_calc and other_peak.is_peak==False:
                    if one_peak.to_check:
                        print_raport(f"{residue} - This residue doesn't have reference peak and check peak position!!!")
                    else: print_raport(f"{residue} - This residue doesn't have reference peak")
                elif one_peak.is_gamma_calc==False and other_peak.is_peak:
                    if one_peak.to_check:
                        print_raport(f"{residue} - This residue doesn't have reference value of gamma and check peak position!!!")
                    else: print_raport(f"{residue} - This residue doesn't have reference value of gamma")
                elif one_peak.is_overlap:
                    if one_peak.to_check:
                        print_raport(f"{residue} - This residue overlaps with: {one_peak.overlap_peaks} and check peak position!!!")
                    else: print_raport(f"{residue} - This residue overlaps with: {one_peak.overlap_peaks}")
                else:
                    # residue = str(other_peak.aa_number)+other_peak.aa_name
                    if one_peak.to_check:
                        print_raport(f"{residue} - This residue doesn't have reference value of gamma or reference peak and check peak position!!!")
                    else: print_raport(f"{residue} - This residue doesn't have reference value of gamma or reference peak")
    
    
    def WriteCCRRate_all_info_new(self, file_type:str):
        """ Methods for writing info for one CCR rate.

        file_type: "csv" or "txt"
        """
        add = self.Additional_text()
        max_lenth_intens = 0
        blank_text = '--'
        max_lenth_dict = {}

        for p in self._peaks:
            if p.is_peak:
                if len(str(p.peak_intens[0]))>max_lenth_intens:
                    max_lenth_intens=len(str(p.peak_intens[0]))

        if file_type == 'csv':
            new_list = "{}{}{}_CCRrates.csv".format(RaportDir,self._CCR_name,add)
        elif file_type == 'txt':
            new_list = "{}{}{}_CCRrates.list".format(RaportDir,self._CCR_name,add)

        with open(new_list, mode='w', newline='') as new_file:
            # headers = ['AA','peak position','intensity in reference','intensity in transfer','CCR rate','Comments']
            headers = ['AA']
            W_list = []
            for i in range(1,self._n_dim+1):
                headers.append(deepcopy('w{}'.format(i)))
                W_list.append(deepcopy('w{}'.format(i+1)))
            headers.extend(deepcopy(['intensity in reference',"Peak uncertainty (reference)",'intensity in transfer',"Peak uncertainty (transfer)",'CCR rate','Uncertainty','Reference Gamma', 'Ix Theor', 'Ix/Ia', 'Ix/Ia without NS','Comments']))

            # max_lenth_discrip = 0
            # for sentence in headers:
            #     if len(sentence)>max_lenth_discrip:
            #         max_lenth_discrip=len(sentence)

            for sentence in headers:
                if sentence in W_list:
                    max_lenth_dict[sentence] = 9
                else:
                    max_lenth_dict[sentence] = len(sentence)+2

            if file_type == 'csv':
                writer_row = csv.writer(new_file,delimiter=",")
                writer_row.writerow(['reference list name =',self._ref_name[0]])
                writer_row.writerow(['transfer list name =',self._trans_name[0]])
                writer_row.writerow(['Tc =', self._tc_vol])
                writer_row.writerow(['$R^2$ =  =', self.r2])
                writer_row.writerow(['---','---','---','---','---','---','---','---','---'])
                # writer_row.writerow(['AA','{:^{sentence_len1}}'.format('peak position',sentence_len1=7*s_dim),'intensity in reference','intensity in transfer','CCR rate','Comments'])
                writer = csv.DictWriter(new_file, fieldnames=headers, delimiter=",")
                writer.writeheader()
            elif file_type == 'txt':
                print ("\tauto list name = {}\tcross list name = {}\tTc = {}\t$R^2$ = {}".format(self._ref_name[0], 
                                                                                     self._trans_name[0], 
                                                                                     self._tc_vol,
                                                                                     self.r2), 
                                                                                     file=new_file) 
                for txt_i in headers:
                    print("{:^{sentence_len}}".format(txt_i,sentence_len=max_lenth_dict[txt_i]), end="\t", file=new_file)
                print("",file=new_file)

            for one_peak in self._peaks:
                one_row = {}
                peak_descr = str(one_peak.aa_number)+one_peak.aa_name
                one_row["AA"] = peak_descr
                if one_peak.is_peak:
                    for i in range(self._n_dim):
                        one_row['w{}'.format(i+1)] = '{:.3f}'.format(one_peak.peak_pos[0][i])
                    one_row["intensity in reference"] = "{:.2e}".format(one_peak.peak_intens[0])
                    one_row["intensity in transfer"] = "{:.2e}".format(one_peak.peak_intens[1])
                    one_row["Peak uncertainty (reference)"] = "{:.2e}".format(one_peak.peak_uncertainty[0])
                    one_row["Peak uncertainty (transfer)"] = "{:.2e}".format(one_peak.peak_uncertainty[1])
                if self._CCR_name == "CCR_10":
                    pass
                if one_peak.is_ccr_rate and not one_peak.is_overlap and not one_peak.ccrrate_calculation_error:
                    one_row["CCR rate"] = '{:.4f}'.format(one_peak.ccr_rate)
                else:
                    one_row["CCR rate"] = blank_text

                if one_peak.ccrrate_error_value != -99999.0:
                    one_row['Uncertainty'] = '{:.4f}'.format(one_peak.ccrrate_error_value)            # ERROR GAMMA INFO, 9.11.2023
                else:
                    one_row['Uncertainty'] = blank_text
                
                one_row["Comments"] = ""

                if one_peak.aa_name == "G":
                    one_row["Comments"] = " | Glycine - no CCR rate"
                if one_peak.to_check:
                    one_row["Comments"] += " | Check peak position!"
                if one_peak.is_overlap == True:
                    one_row['Comments'] += " | This peak overlap with: {}".format(one_peak.overlap_peaks)
                # if one_peak.is_peak == False and one_peak.is_ccr_rate:
                #     one_row['Comments'] = "CCR rate value from other peak:"
                if one_peak.ccrrate_calculation_error:
                    one_row["Comments"] += ' | Peak intens error: Ix/Ia = {:.4f}'.format(one_peak.ccr_rate)
                if one_peak.is_peak == False:
                    one_row["Comments"] += ' | No visible peak'
                if one_peak.fatal_error:
                    one_row["Comments"] += ' | Fatal Error - experimental data faw away from reference'
                    

                if one_peak.is_gamma_calc and one_peak.is_ccr_rate:
                    one_row['Reference Gamma'] = "{:.4f}".format(one_peak.gamma_ref)
                else:
                    one_row['Reference Gamma'] = blank_text

                if one_peak.is_theor_values:
                    one_row['Ix Theor'] = "{:.2e}".format(one_peak.Ix_theor)
                    one_row['Ix/Ia'] = "{:.2f}".format(one_peak.Ix_theor_Ia_ratio)
                    one_row['Ix/Ia without NS'] = "{:.2f}".format(one_peak.Ix_theor_Ia_ratio_without_NS)
                else:
                    one_row['Ix Theor'] = blank_text
                    one_row['Ix/Ia'] = blank_text
                    one_row['Ix/Ia without NS'] = blank_text
                
                if file_type == 'csv':
                    writer.writerow(one_row)
                elif file_type == 'txt':
                    for txt_i in headers:
                        if txt_i in one_row:
                            print("{:^{sentence_len}}".format(one_row[txt_i],sentence_len=max_lenth_dict[txt_i]), end="\t", file=new_file)
                        else:
                            print("{:^{sentence_len}}".format(blank_text,sentence_len=max_lenth_dict[txt_i]), end="\t", file=new_file)

                    print("",file=new_file)
        # print ("All info from peaks for {} rates are in file: {}".format(ccr_name, new_list), file=RaportBox)
        RaportBox.write("All info from peaks for {} rates are in file: {}\n".format(self._CCR_name, new_list))
        return


class CCR_SymRec(CCRClass):
    
    def __init__(self,exp_dict):
        super().__init__(exp_dict)   

    def Read_peaklist(self, peaklist_directory, points_mode=False):
        NameFlag = [False, False, False, False]
        if points_mode == False:
            list_of_names_ends = [".list","_new_ppm.list"]
        else:
            list_of_names_ends = ["_points.list","_new_points.list"]
        peak_list_basic_names = [self._ref_name[0], self._ref_name[1], 
                                 self._trans_name[0], self._trans_name[1]]
        peak_list_names = ["None","None","None","None"]
        for indexlv, list_verson in enumerate(peak_list_basic_names):
            for i, l in enumerate(list_of_names_ends):
                peaklistfile = f"{peaklist_directory}/peak_lists/{list_verson}{l}"
                try:
                    f=open(peaklistfile)
                    NameFlag[indexlv] = True
                    peak_list_names[indexlv] = peaklistfile
                    f.close
                except FileNotFoundError:
                    RaportBox.write("There is no such file or directory: {}\n".format(peaklistfile))
        
        if NameFlag[0] == True and NameFlag[1] == True:
            RaportBox.write("\n\nLIST:{}\n".format(peak_list_names))
            pl_a_1 = open(peak_list_names[0], 'r')
            pl_a_2 = open(peak_list_names[1], 'r')
            pl_x_1 = open(peak_list_names[2], 'r')
            pl_x_2 = open(peak_list_names[3], 'r')
            self._is_peaklist = True
            p_lines_a_1 = pl_a_1.readlines()
            p_lines_a_2 = pl_a_2.readlines()
            p_lines_x_1 = pl_x_1.readlines()
            p_lines_x_2 = pl_x_2.readlines()
            if points_mode:
            # if points_mode and "Average noise level" in p_lines_a[1] and p_lines_x[1]:
                try:
                    self._noise = [float(p_lines_a_1[1].split("=")[1]), float(p_lines_a_2[1].split("=")[1]), 
                                float(p_lines_x_1[1].split("=")[1]), float(p_lines_x_2[1].split("=")[1])]
                    print_raport(f"auto_noise = {self._noise[0]} and {self._noise[1]}, cross_noise = {self._noise[2]} and {self._noise[3]}")
                except:
                    if len(self._noise)>0:
                        print_raport(f"noise from input file: auto_noise = {self._noise[0]} and {self._noise[1]}, cross_noise = {self._noise[2]} and {self._noise[3]}")
                    else:
                        print_raport(f"No info about noise level")
                # RaportBox.write(f"auto_noise = {self._noise[0]}, cross_noise = {self._noise[1]}")
            for indexl, line in enumerate(p_lines_a_1):
                if indexl > 1 and len(line) > 1:
                    item_a_1 = line.split()
                    item_a_2 = p_lines_a_2[indexl].split()
                    item_x_1 = p_lines_x_1[indexl].split()
                    item_x_2 = p_lines_x_2[indexl].split()
                    # Reading description
                    if item_a_1[0] != item_x_1[0] or item_a_2[0] != item_x_2[0] or  item_a_1[0] != item_x_2[0] or item_a_2[0] != item_x_1[0]:
                        print ("Error: reference ({} / {}) and transfer ({} / {}) peak lists are not compatible:\n{}/{} ? {}/{}".format(peak_list_names[0], peak_list_names[1],
                                                                                                                             peak_list_names[2], peak_list_names[3], 
                                                                                                                            item_a_1[0], item_a_2[0],
                                                                                                                            item_x_1[0], item_x_2[0]))
                        sys.exit()
                    description = item_a_1[0]
                    aminoacids_number = self.read_aminoacid_number(description)
                    RaportBox.write(f"aminoacids_number: {aminoacids_number,description}\n")
                    aminoacids_number = int(aminoacids_number)
                    for res in self._peaks: 
                        if aminoacids_number == res.aa_number:
                            res.is_peak = True
                            res.descript = deepcopy(description)
                            if points_mode == True:  # for points mode
                                for i in range(1,self._n_dim+1):
                                    res.peak_pos_points[0].append(deepcopy(int(item_a_1[i])))
                                    res.peak_pos_points[1].append(deepcopy(int(item_a_2[i])))
                                    res.peak_pos_points[2].append(deepcopy(int(item_x_1[i])))
                                    res.peak_pos_points[3].append(deepcopy(int(item_x_2[i])))
                                    # print ("residues append --->", residue.peak_pos_points)
                            else: 
                                for i in range(1,self._n_dim+1):
                                    res.peak_pos[0].append(deepcopy(float(item_a_1[i])))
                                    res.peak_pos[1].append(deepcopy(float(item_a_2[i])))
                                    res.peak_pos[2].append(deepcopy(float(item_x_1[i])))
                                    res.peak_pos[3].append(deepcopy(float(item_x_2[i])))
                    # Reading peak intensity
                            if len(item_a_1)>self._n_dim+1:
                                RaportBox.write("reference = {} and {}, transfer = {} and {}\n".format(item_a_1[self._n_dim+1],item_a_2[self._n_dim+1],
                                                                                               item_x_1[self._n_dim+1], item_x_2[self._n_dim+1]))
                                res.peak_intens = [ float(item_a_1[self._n_dim+1]) , float(item_a_2[self._n_dim+1]),
                                                   float(item_x_1[self._n_dim+1]) , float(item_x_2[self._n_dim+1]) ]
        else:
            print_raport ("Missing file for {} / {} and {} / {}! {}\n".format(self._ref_name[0], self._ref_name[1], 
                                                                              self._trans_name[0], self._trans_name[1], 
                                                                              peak_list_names))
            self._is_peaklist = False
       
    def calc_uncertainty_value(self,peak_number):
        Ia = self._peaks[peak_number].peak_intens[0]
        Ix = self._peaks[peak_number].peak_intens[1]
        NSa = self._ns[0]
        NSx = self._ns[1]
        noiseIa = self._noise[0]
        noiseIx = self._noise[1]
        uncertainty_Ia = self._peaks[peak_number].peak_uncertainty[0]
        uncertainty_Ix = self._peaks[peak_number].peak_uncertainty[1]
        print_raport(f"{self._peaks[peak_number].descript}: Peak_uncertainty for CalcErrorValue: {uncertainty_Ia}, {uncertainty_Ix}")
        
        squareDerivativeArcTanh = math.pow(1/(1-math.pow((Ix*NSa)/(Ia*NSx),2)),2)
        sumOfSquaresOfPeakHightDeviation = math.pow(noiseIx/Ia,2) + math.pow(noiseIa*Ix/math.pow(Ia,2),2)
        squareError = math.pow(1/self._tc_vol,2) * squareDerivativeArcTanh * math.pow(NSa/NSx,2) * sumOfSquaresOfPeakHightDeviation
        self._peaks[peak_number+self._CCR_pos].ccrrate_error_value = math.sqrt(squareError)

        print_raport("Calculating gamma uncertainty for peak number: {} \n squareDerivativeArcTanh = {} \n sumOfSquaresOfPeakHightDeviation = {}\n squareError = {}\n \t\t--> {}\n".format(peak_number,
                                                                                                                                                        squareDerivativeArcTanh,
                                                                                                                                                        sumOfSquaresOfPeakHightDeviation,
                                                                                                                                                        squareError,
                                                                                                                                                        math.sqrt(squareError)))
    
    def calc_ccr_rate(self):
        for indexp, peak in enumerate(self._peaks):
            if not peak.is_overlap and peak.is_peak:
                res_num = indexp+self._CCR_pos
                other_peak = self._peaks[res_num]
                self.calc_uncertainty_value(indexp)
                try:
                    ccr_rate_vol = atanh(math.sqrt((peak.peak_intens[2]*self._ns[0]*peak.peak_intens[3]*self._ns[1])/(peak.peak_intens[0]*self._ns[2]*peak.peak_intens[1]*self._ns[3])))/self._tc_vol
                    other_peak.is_ccr_rate = True
                except:
                    ccr_rate_vol = ((peak.peak_intens[2]*self._ns[0]*peak.peak_intens[3]*self._ns[1])/(peak.peak_intens[0]*self._ns[2]*peak.peak_intens[1]*self._ns[3])) #"atanh(x) - x shoudl be beetween -1 and 1"
                    other_peak.ccrrate_calculation_error = True
                other_peak.ccr_rate = ccr_rate_vol
                
        return

    def calc_theor_Ix(self):
        for indexp, one_peak in enumerate(self._peaks):
            res_num = indexp-self._CCR_pos
            if len(self._peaks)> res_num >= 0:
                other_peak = self._peaks[res_num]
                if one_peak.is_gamma_calc and other_peak.is_peak:
                    residue = str(other_peak.aa_number)+other_peak.aa_name
                    print_raport(residue)
                    auto_peak_intensity = other_peak.peak_intens[0]
                    other_peak.Ix_theor = float(math.tanh(one_peak.gamma_ref*self._tc_vol)*auto_peak_intensity*self._ns[1]/self._ns[0])
                    if self._CCR_name == "CCR_5" or self._CCR_name == "CCR_6":
                        other_peak.Ix_theor = -1 * other_peak.Ix_theor
                    other_peak.Ix_theor_Ia_ratio = one_peak.Ix_theor/auto_peak_intensity
                    other_peak.Ix_theor_Ia_ratio_without_NS = one_peak.Ix_theor_Ia_ratio*self._ns[0]/self._ns[1]
                else:
                    residue = str(other_peak.aa_number)+other_peak.aa_name
                    print_raport("This residue ({}) doesn't have reference value of gamma or reference peak".format(residue))


    def read_peak_files(self,peaklist_directory:str,seq_dict:dict):
        self.Prepare_peaklist(seq_dict)
        for one_peak in self._peaks:
            one_peak.peak_uncertainty = [0.0, 0.0, 0.0, 0.0]
            one_peak.peak_pos = [[],[],[],[]]
            one_peak.peak_pos_points = [[],[],[],[]]
        self.Read_peaklist(peaklist_directory, points_mode=True)
        self.Read_peaklist(peaklist_directory, points_mode=False)
        self.Read_peak_uncertainty(peaklist_directory, self._ref_name[0], version=0)
        self.Read_peak_uncertainty(peaklist_directory, self._ref_name[1], version=1)
        self.Read_peak_uncertainty(peaklist_directory, self._trans_name[0], version=2)
        self.Read_peak_uncertainty(peaklist_directory, self._trans_name[1], version=3)




class CResidue:
    def __init__(self,aa_num,aa_name):
        self.aa_name = aa_name
        self.aa_number = aa_num              # the number of the amino acid to which this peak corresponds
        self.is_peak = False       # information about the presence of a peak

        self.peak_pos = []            # chemical shifts for all nuclei of peak, length depends on dimentionality
        self.peak_pos_points = []
        self.peak_intens = []      # peak hight in reference [0] and transfer [1] version
        self.descript = ""              # description of peak which is in first column in Sparky-like peak list
        # self.peak_seq_pos_a = []
        self.is_overlap = False       # information about the presents of overlapping of a peak: maybe, yes, no   
        self.overlap_peaks = {}       #     
        self.peak_uncertainty = [] 
        self.to_check = False

        self.is_ccr_rate = False
        self.ccr_rate = -99999.0
        self.ccrrate_error_value = -99999.0
        self.ccrrate_calculation_error = False

        self.is_gamma_calc = False      # information about the presence of a reference gamma
        self.gamma_ref = -99999.0       # value of gamma (CCR rate) which will use as reference data for experimental data, value additional file with calculated previosly gammas
        self.fatal_error = False        # experimental data faw away from reference (more than 20% range of rates)
        self.is_theor_values = False
        self.Ix_theor = 'nan'
        self.Ix_theor_Ia_ratio = 'nan'
        self.Ix_theor_Ia_ratio_without_NS = 'nan'


    # def check_if_fatal_error(self):
    #     if self.is_gamma_calc:
    #         dist = math.sqrt((self.ccr_rate - self.gamma_ref )**2)
    #         if dist > 5:
    #             self.fatal_error = True
            # quotient = self.ccr_rate/self.gamma_ref
            # if 1.25 < quotient < 0.75:
            #     self.fatal_error = True

"""Small functions"""


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



"""Other functions"""

                

def ReadSequnce(seq_file_name:Path) -> dict:
    seq_dict = {}
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
                        aa_no += 1
                        seq_dict[deepcopy(aa_no)] = deepcopy(aa)
            print_raport ("FASTA sequence: {}".format(FASTAseq))
    return seq_dict

def Read_Reference_Gamma(gamma_cal_file_name:Path)->dict[dict[str:float]]:
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
                # print ("gamma_file_dict: AA - {}, NUM - {}".format(gamma_file_dict["AA"][i],gamma_file_dict["AA"][i][:-3]))
                gamma_file_dict["seq_num"].append(deepcopy(gamma_file_dict["AA"][i][:-3]))
    except:
        if gamma_file_dict["AA"][0][:-1].isnumeric():
            gamma_file_dict["seq_num"]=[gamma_file_dict["AA"][0][:-1]]
            for i in range(1,len(gamma_file_dict["AA"])):
                gamma_file_dict["seq_num"].append(deepcopy(gamma_file_dict["AA"][i][:-1]))
    # print(gamma_file_dict)
    return gamma_file_dict


def LRegression_expresion(x,y):
    slope, intercept, r, p, std_err = stats.linregress(x, y)
    matching_factor_for_y_x = (1-slope)/std_err #if it is <1 is good, >1 means that we have weak corellation x=y 
    r2 = r**2
    lr_exp = '{:.3f}*x + {:.3f}'.format(slope,intercept)
    RaportBox.write(f"\nlinear regression: {lr_exp}, r2: {r**2:.3f}\n p: {p}, std_err:{std_err}\n")
    lr_exp_y_list = []
    for i in x:
        lr_exp_y_list.append(deepcopy(slope*i+intercept))
    return lr_exp, r2, slope, intercept, matching_factor_for_y_x

def WeightedLRegression_expresion(x,y,uncertainty_val):
    from sklearn.linear_model import LinearRegression
    regr = LinearRegression()
    # X = np.array([[i,i] for i in x if x!='nan'])
    # Y = np.array([[i,y[index]] for index, i in enumerate(x)])
    X = np.array([[i] for i in x if x!='nan'])
    Y = np.array([[i] for i in y if x!='nan'])
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
    RaportBox.write(f"weighted linear regression: {lr_exp}, r2: {r2:.3f}")
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
    RaportBox.write("SSR = {:.0f}; SSE = {:.0f}; SST = {:.0f}; SSR+SSE = {:.0f}; R2_1 = {:.4f}; R2_2 = {:.4f}".format(SSR, SSE, SST, SSR+SSE, R2_val, R2_val2))
    return R2_val


def WeightedLRegression_expresion_by_hand(x:list,y:list,uncertainty_val:list)->dict["equation":str, 
                                                                                    "slope":float, 
                                                                                    "slope_uncertainty":float, 
                                                                                    "intercept":float, 
                                                                                    "intercept_uncertainty":float, 
                                                                                    "factor_a_for_y_x":float, 
                                                                                    "factor_b_for_y_x":float, 
                                                                                    "r2":float]:
    """
    output: dict["equation":str, 
                            "slope":float, 
                            "slope_uncertainty":float, 
                            "intercept":float, 
                            "intercept_uncertainty":float, 
                            "factor_a_for_y_x":float, 
                            "factor_b_for_y_x":float, 
                            "r2":float]
    """
    if len(x)>1:
        RaportBox.write("\nby Hand\n")
        list_of_invert_squer_uncertainty = []
        list_of_avg_x = []
        list_of_avg_square_x = []
        list_of_avg_y = []
        list_of_x_avgx_y_avgy_div_uncertainty = []
        list_of_x_avgx_squer_div_uncertainty = []

        for i in range(len(x)):
            if str(x[i])!='nan' and str(y[i])!='nan' and str(uncertainty_val[i])!='nan':
                list_of_invert_squer_uncertainty.append(deepcopy(1/(uncertainty_val[i]**2)))
                list_of_avg_x.append(deepcopy(x[i]/(uncertainty_val[i]**2)))
                list_of_avg_y.append(deepcopy(y[i]/(uncertainty_val[i]**2)))
                list_of_avg_square_x.append(deepcopy((x[i]**2)/(uncertainty_val[i]**2)))

        sum_of_invert_squer_uncertainty = sum(list_of_invert_squer_uncertainty)
        avg_x = sum(list_of_avg_x)/sum_of_invert_squer_uncertainty
        avg_y = sum(list_of_avg_y)/sum_of_invert_squer_uncertainty
        avg_square_x = sum(list_of_avg_square_x)/sum_of_invert_squer_uncertainty
        flag = False
        for i in range(len(x)):
            if str(x[i])!='nan' and str(y[i])!='nan' and str(uncertainty_val[i])!='nan':
                flag = True
                list_of_x_avgx_y_avgy_div_uncertainty.append(deepcopy((x[i]-avg_x)*(y[i]-avg_y)/(uncertainty_val[i]**2)))
                list_of_x_avgx_squer_div_uncertainty.append(deepcopy(((x[i]-avg_x)/uncertainty_val[i])**2))

        sum_of_x_avgx_y_avgy_div_uncertainty = sum(list_of_x_avgx_y_avgy_div_uncertainty)
        sum_of_x_avgx_squer_div_uncertainty = sum(list_of_x_avgx_squer_div_uncertainty)
        # print_raport(f"sum_of_x_avgx_y_avgy_div_uncertainty: {sum_of_x_avgx_y_avgy_div_uncertainty}, sum_of_x_avgx_squer_div_uncertainty: {sum_of_x_avgx_squer_div_uncertainty}")
        slope = sum_of_x_avgx_y_avgy_div_uncertainty/sum_of_x_avgx_squer_div_uncertainty  
        intercept = avg_y-(slope*avg_x)
        slope_uncertainty = math.sqrt(1/(sum_of_x_avgx_squer_div_uncertainty))
        intercept_uncertainty = math.sqrt((avg_square_x)/sum_of_x_avgx_squer_div_uncertainty)

        # print_raport(f"slope: {slope:.3f},intercept:{intercept:.3f},slope_uncertainty:{slope_uncertainty:.5f},intercept_uncertainty:{intercept_uncertainty:.5f}")
        # lr_exp = '{:.2f}*x + {:.2f}'.format(slope,intercept)
        # # print_raport(f"weighted linear regretion: {lr_exp}")
        # factor_a_for_y_x = (1-abs(slope))/slope_uncertainty
        # factor_b_for_y_x = (1-abs(intercept))/intercept_uncertainty
        # print_raport(f"1-a/delta_a = {factor_a_for_y_x}")
        
        y_from_line = []
        for one_x in x:
            if str(one_x)!='nan':
                y_from_line.append(deepcopy(one_x*slope+intercept))
        # r2 = calc_R2(y,y_from_line)
        
        # return_dict = {"equation":'{:.2f}*x + {:.2f}'.format(slope,intercept), 
        #                "slope":slope, 
        #                "slope_uncertainty":slope_uncertainty, 
        #                "intercept":intercept, 
        #                "intercept_uncertainty":intercept_uncertainty, 
        #                "factor_a_for_y_x":(1-abs(slope))/slope_uncertainty, 
        #                "factor_b_for_y_x":(1-abs(intercept))/intercept_uncertainty, 
        #                "r2":calc_R2(y,y_from_line)}
        return {"equation":'{:.2f}*x + {:.2f}'.format(slope,intercept), 
                "slope":slope, 
                "slope_uncertainty":slope_uncertainty, 
                "intercept":intercept, 
                "intercept_uncertainty":intercept_uncertainty, 
                "factor_a":(1-abs(slope))/slope_uncertainty, 
                "factor_b":(1-abs(intercept))/intercept_uncertainty, 
                "r2":calc_R2(y,y_from_line)}
    else:
        print ("Too small number of point to calculate Regression - minimal is 2")

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

def set_nrow_ncol(number_of_experiment):
    if 13<=number_of_experiment<=16:
        nrows=4
        ncols=4
    elif 10<=number_of_experiment<=12:
        nrows=3
        ncols=4
    elif number_of_experiment==9:
        nrows=3
        ncols=3
    elif number_of_experiment==8 or number_of_experiment==7:
        nrows=2
        ncols=4
    elif number_of_experiment==6 or number_of_experiment==5:
        nrows=2
        ncols=3
    elif number_of_experiment==4 or number_of_experiment==3:
        nrows=2
        ncols=2
    else:
        nrows=3
        ncols=3
    return nrows,ncols
"""                      MAIN PROGRAM                      """

if __name__ == "__main__":
    RaportBox = open(RaportBoxFile,"w")
    print_raport("\n=== Reading input files ===")
    SequenceDict = ReadSequnce(seq_file_name)

    if refgammaFlag:
        gamma_ref_dict = Read_Reference_Gamma(gamma_cal_file_name)
    else: gamma_ref_dict = {}

    ExperimentsSet = CCRSet(file_directory=file_directory, 
                            peaklist_directory=peaklist_directory,
                            exp_file=exp_file_name,
                            seq_dict=SequenceDict,
                            ref_flag=refgammaFlag)
    ExperimentsSet.CalcCCRrates()
    plot_together_list = []
    
    if refgammaFlag:
        ExperimentsSet.Compere_with_reference(ref_dict=gamma_ref_dict)

    for key in ExperimentsSet.to_compere_dict:
        exp = ExperimentsSet.ccr_set[ExperimentsSet.to_compere_dict[key][-1]]
        plot_together_list.append(deepcopy(exp))
    
    if len(plot_together_list)>0 and refgammaFlag:
        ExperimentsSet.plot_gamma_gamma_all_together(plot_together_list, style=output_style)
    
    ExperimentsSet.Compere_diff_ver_exp()

    RaportBox.close()


    # TO DO
    # dodac ograniczenie dla numerów reszt dla których liczymy 
    # dodac dodawanie daty i adnotacji do wykresow 
