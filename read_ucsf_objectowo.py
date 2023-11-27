#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

"""
Created on Feb 15 12:55 2022

@author: Paulina Bartosinska-Marzec
"""
""" The script for reading spectra in ucsf format and outputting values for given points in the spectrum """


from copy import deepcopy
from random import randint
import struct
import sys
from itertools import product
import os
import numpy
import math
import time
import argparse #https://docs.python.org/3/library/argparse.html
from pathlib import Path
from BioNMRdictionary import *
import logging


parser = argparse.ArgumentParser(
                    prog='read_ucsf',
                    description="read_ucsf script after reading ucsf file and peak list (in Sparky format), check peak intensity of peak and if it is not in the highest position - move it. \n After this script prints peak lists in ppm value and points value.",
                    epilog='Text at the bottom of help')

parser.add_argument("filename", metavar="ucsf_path", type=Path, 
                    help="path to UCSF file")

parser.add_argument("peak_list", metavar="peak_list_path", type=Path, 
                    help="path to peak list in SPARKY format")

parser.add_argument("-np", "--npoints", dest='Number_of_points_for_noise', type=int, default=10, 
                    help="to change number of points for calculate noise level: N^(spectra dimentionality + 1); normally is N = 10")

parser.add_argument("-pl", "--peaklevel", dest='peak_level', type=float, default=0.0,
                    help="if you know level when starting appear, add this with scientific numer notation e.g. 1e+7")

parser.add_argument("-nrm", "--noRemove", dest='noRemoveFlag', action="store_true", default=False,
                    help="add this if you do not want remove invisible peaks")

parser.add_argument("-op", "--onlypoints", dest='OnlyPoints', action="store_true",  default=False,
                    help="add this if you want only change ppm value to points value")

parser.add_argument("-o", "--output_name", type=Path, 
                    help="add this if you want specific output name")

parser.add_argument("-n", "--noise", dest='OnlyNoise', action="store_true",  default=False,
                    help="add this if you want calculate only noise level")

args = parser.parse_args()



if not os.path.exists(args.filename): 
    print("There is not such file like this: ", args.filename)
    sys.exit(1)
if not os.path.exists(args.peak_list): 
    print("There is not such file like this: ", args.peak_list)
    sys.exit(1)


if args.peak_level:
    UserPeakLevelFlag = True
else: 
    UserPeakLevelFlag = False

if args.noRemoveFlag:
    print ("Peak centering without removing invisible peaks")
if args.OnlyPoints:
    print ("OnlyPoints is on")    
if args.OnlyNoise:
    print ("Only calculate noise level option is on")    

if args.output_name:
    print(args.output_name)
    peak_list_name = os.path.basename(args.output_name)
    if "." in str(peak_list_name):
        point_index = str(peak_list_name).index(".")
        print("index . : {}".format(point_index))
        peak_list_name = peak_list_name[:point_index]

    peak_list_dir_new = os.path.dirname(args.output_name)+peak_list_name+"_list"
    if not os.path.exists(peak_list_dir_new):
        os.mkdir(peak_list_dir_new)
    print ("Output file path: {}/{}.list".format(peak_list_dir_new,peak_list_name))
else: 
    peak_list_name = os.path.basename(args.peak_list)
    if "." in str(peak_list_name):
        point_index = str(peak_list_name).index(".")
        peak_list_name = peak_list_name[:point_index]

    peak_list_dir_new = os.path.dirname(args.peak_list)+peak_list_name+"_list"
    if not os.path.exists(peak_list_dir_new):
        os.mkdir(peak_list_dir_new)
    print ("Output file path: {}/{}".format(peak_list_dir_new,peak_list_name))




def print_raport(text_to_write,in_terminal=True):
    with open('{}/info.txt'.format(peak_list_dir_new), 'a') as txtfile:
        if in_terminal:
            print (text_to_write)
        txtfile.write(text_to_write+"\n")

class CSpectrum:
    def __init__(self,filename,peak_list):
        self.nucl_name = []
        self.points_num = []
        self.tile_size = []
        self.n_tiles=[]      # number of tiles along axis [direct, 1, 2, 3] ,# point number with peak in tile  
        self.spectr_fq =[]
        self.sw = []
        self.sw_ppm = []
        self.data_center = []
        self.spectra_dimentionality = 0
        self.peak_level = 0.0
        self.noise_level = 0.0
        
        self.ReadSpecraParamiters(filename)
        self.ChangeHz2ppm()
        
        self.peaks = self.read_peaklist(peak_list)
        
    def ReadSpecraParamiters(self,filename):
        with open(filename, "rb") as ucsf_file:
            """Read spectra dimentionality"""
            ucsf_file.seek(10)
            ucsf_data = ucsf_file.read(1)                                   # 1 byte per information 
            self.spectra_dim = int.from_bytes(ucsf_data, byteorder='big')         # convert bytes to integer value 
            print ("Spectra dimensiolity:", self.spectra_dim)

            """Read spectra format"""
            ucsf_file.seek(13)
            ucsf_data = ucsf_file.read(1)                                   # 1 byte per information 
            spectra_format = int.from_bytes(ucsf_data, byteorder='big')         # convert bytes to integer value  

            """Read nucleus name"""
            for i in range(0,self.spectra_dim):
                ucsf_file.seek(180+128*i)                                       # nucleus name (1H, 13C, 15N, 31P); first is after 180 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(6)                                   # 6 byte per information 
                text_data = ucsf_data.decode('utf-8')                           # convert bytes to string    
                # print ("nucleus name", ucsf_data, "\t", text_data)
                # print (text_data)
                self.nucl_name.append(deepcopy(text_data))
                # print (self.nucl_name[0])
            self.nucl_name.insert(0,self.nucl_name.pop())
            # print ("tutaj", self.nucl_name[0], self.nucl_name[1],self.nucl_name[2],self.nucl_name[3])
            # for ii in range(len(self.nucl_name)):
            #     print("bbb", self.nucl_name[ii])
            print ("Nucleus names:",*self.nucl_name)

            """Read number of points per asix"""
            for i in range(0,self.spectra_dim):
                ucsf_file.seek(188+128*i)                                       # integer number of data points along axis; first is after 188 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(4)                                   # 4 byte per information 
                bytes2int = int.from_bytes(ucsf_data, byteorder='big')          # convert bytes to integer value  
                # print ("self.points_num", ucsf_data, "\t", bytes2int)
                self.points_num.append(deepcopy(bytes2int)) 
            self.points_num.insert(0,self.points_num.pop())
            print ("Number of points in axis:",*self.points_num)

            """Read integer tile size along this axis"""
            for i in range(0,self.spectra_dim):
                ucsf_file.seek(196+128*i)                                       # integer tile size along this axis; first is after 196 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(4)                                   # 4 byte per information 
                bytes2int = int.from_bytes(ucsf_data, byteorder='big')          # convert bytes to integer value  
                # print ("tile size", ucsf_data, "\t", bytes2int)
                self.tile_size.append(deepcopy(bytes2int))
            self.tile_size.insert(0,self.tile_size.pop())
            # print (self.tile_size)

            """Read float spectrometer frequency for this nucleus (MHz) """
            for i in range(0,self.spectra_dim):
                ucsf_file.seek(200+128*i)                                       # float spectrometer frequency for this nucleus (MHz) ; first is after 196 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(4)                                   # 4 byte per information 
                [bytes2float] = struct.unpack('>f', ucsf_data)                  # convert bytes to float value
                self.spectr_fq.append(deepcopy(bytes2float))
            self.spectr_fq.insert(0,self.spectr_fq.pop())
            print ("Spectrometer frequency (MHz): ",' '.join("{:.2f}".format(x) for x in self.spectr_fq))

            """Read float spectral width (Hz)"""
            for i in range(0,self.spectra_dim):
                ucsf_file.seek(204+128*i)                                       # float spectral width (Hz) ; first is after 196 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(4)                                   # 4 byte per information 
                [bytes2float] = struct.unpack('>f', ucsf_data)                  # convert bytes to float value
                self.sw.append(deepcopy(bytes2float))
            self.sw.insert(0,self.sw.pop())
            print ("Spectral width (Hz):", ' '.join("{:.0f}".format(x) for x in self.sw))


            """Read float center of data (ppm)"""
            for i in range(0,self.spectra_dim):
                ucsf_file.seek(208+128*i)                                       # float center of data (ppm) ; first is after 196 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(4)                                   # 4 byte per information 
                [bytes2float] = struct.unpack('>f', ucsf_data)                  # convert bytes to float value
                self.data_center.append(deepcopy(bytes2float))
            self.data_center.insert(0,self.data_center.pop())
            # print ("self.data_center",self.data_center)
    
    def ChangeHz2ppm(self):
        """Change spectral width from Hz to ppm"""
        for p in range(len(self.sw)):
            sw_sparky = self.sw[p]-(self.sw[p]/self.points_num[p])
            self.sw_ppm.append(deepcopy(sw_sparky/self.spectr_fq[p]))
        print ("Spectral width (ppm):", ' '.join("{:.1f}".format(x) for x in self.sw_ppm))
        # print ("self.spectr_fq", self.spectr_fq)
        # print ("self.data_center",self.data_center)
        # print ("self.sw =", self.sw, "\nsw_ppm =", self.sw_ppm)
        
        for i in range (0,self.spectra_dim):
            if self.points_num[i]%self.tile_size[i]==0:
                self.n_tiles.append(deepcopy(int(self.points_num[i]/self.tile_size[i])))
            else: 
                self.n_tiles.append(deepcopy(int(self.points_num[i]/self.tile_size[i])+1))

    def read_peaklist(self, peak_list):
        with open(peak_list, 'r') as pl:
            p_lines = pl.readlines()
            p_list = []
            for indexl, line in enumerate(p_lines):
                if indexl > 1 :
                    p_pos = CPeak(line,self.spectra_dimentionality)
                    p_list.append(deepcopy(p_pos))
                    # print (p_pos.peak_ppm_pos)       
        return p_list
        
    def SetUpPeakLevel(self):
        if UserPeakLevelFlag:
            print_raport("Peak level starts from = {:.2e}".format(args.peak_level))
            self.peak_level = args.peak_level
            self.noise_level = deepcopy(self.peak_level/50)
                
    def CalcNoise(self):
        list_of_noise_points = find_noise_point(self.peaks, self.spectra_dimentionality, self.points_num, args.Number_of_points_for_noise)
        Points_intens = read_intens(list_of_noise_points, args.filename, self.spectra_dimentionality, self.tile_size, self.n_tiles)
        print ("list of points of noise is ready: {} points".format(len(Points_intens)))
        average_noise_level = round(sum(Points_intens)/len(Points_intens), 1)
        sum_of_squares=0.0
        for val in Points_intens:
            sum_of_squares+=(val-average_noise_level)**2
        self.noise_level = math.sqrt(sum_of_squares/len(Points_intens))
        self.peak_level = self.noise_level * 50
    
    


class CPeak():
    def __init__(self,input_line,s_dim):
        self.peak_ppm_pos = []            # chemical shifts for all nuclei of peak, length depends on dimentionality
        self.peak_points_pos = []      # peak position in points
        self.is_visible = True         # if peak height is 50 time biger than noise level = True
        self.peak_intens = 0            # peak height  
        self.descript = ""              
        self.is_center = False          # if peak is in the highest point = True
        self.was_moved = False          # if peak was moved = True
        self.new_points_pos = []         # if peak was moved there will be new position in points
        self.new_ppm_pos = []           # if peak was moved there will be new position in ppm

        if type(input_line) == str:
            input_line = input_line.split()
        if type(input_line) == list:
            self.descript = input_line[0]
            # print ("p_pos.peak_ppm_pos", p_pos.peak_ppm_pos)
            self.peak_ppm_pos.append(deepcopy(float(input_line[s_dim])))
            for i in range(1,s_dim):
                self.peak_ppm_pos.append(deepcopy(float(input_line[i])))
        else:
            pass
            # raise ERROR 

    def calc_peak_points(self, spectrum:CSpectrum, num_format="Integer"):
        for p in range(len(self.peak_ppm_pos)):
            one_dim = ppm2points(spectrum, self.peak_ppm_pos[p],p,num_format)
            self.peak_points_pos.append(deepcopy(one_dim))




"""Calculating functions"""





def Print_Peak_List_ppm(peaks, s_dim, points_num, sw_ppm, data_center):
    new_peak_list = "{0}/{1}_new_ppm.list".format(peak_list_dir_new,peak_list_name)
    print ("new_peak_list", new_peak_list)
    max_lenth_discrip = 0
    for p in peaks:
        if len(p.descript)>max_lenth_discrip:
            max_lenth_discrip=len(p.descript)
    with open(new_peak_list, 'w') as listfile:
        print ("\tAssignment", end="", file=listfile) 
        for i in range(s_dim):
            print ("\tw{}".format(i+1), end="", file=listfile)
        print ("\tData Height\t\t\t\tCentering peak list prepare from{}\n".format(peak_list_name), file=listfile)
        for indexpeak, one_peak in enumerate(peaks):
            if one_peak.is_visible == True:
                print ("{:{sentence_len}}".format(one_peak.descript, sentence_len=max_lenth_discrip), end="\t", file=listfile)
                if one_peak.was_moved == True and one_peak.is_center != False:
                    for i in range(1,s_dim):
                        ppm_position = points2ppm(one_peak.new_points_pos[i], points_num[i],data_center[i],sw_ppm[i])
                        print ("{:.3f}".format(ppm_position), end="\t", file=listfile)
                    ppm_position = points2ppm(one_peak.new_points_pos[0], points_num[0],data_center[0],sw_ppm[0])
                    print ("{:.3f}".format(ppm_position), end="\t", file=listfile)
                        
                else: 
                    for i in range(1,s_dim):
                        print ("{:.3f}".format(one_peak.peak_ppm_pos[i]), end="\t", file=listfile)
                    print ("{:.3f}".format(one_peak.peak_ppm_pos[0]), end="\t", file=listfile)
                print (one_peak.peak_intens, end=" ", file=listfile) 
                
                if one_peak.is_center == False:
                    print ("to check", file=listfile)
                else: print (file=listfile)
    return

def Print_Peak_List_points(peaks, s_dim, type_list):
    new_peak_list = "{}/{}_{}_points.list".format(peak_list_dir_new,peak_list_name,type_list)

    
    with open(new_peak_list, 'w') as listfile:
        print ("\tAssignment", end="", file=listfile) 
        for i in range(s_dim):
            print ("\tw{}".format(i+1), end="", file=listfile)
        if type_list == "orgin":
            print ("\tData Height\t\t\t\t Orgin peaklist in points\n", file=listfile)
            # print ("{} peak list".format(type_list), new_peak_list)
        else:
            print ("\tData Height\t\t\t\tCentering peak list prepare from{}\n".format(peak_list_name), file=listfile)
            print ("{} peak list".format(type_list), new_peak_list)
        for one_peak in peaks:
            if one_peak.is_visible == True:
                print ("{}".format(one_peak.descript), end="\t", file=listfile)
                if one_peak.was_moved == True and one_peak.is_center != False:
                    for i in range(1,s_dim):
                        print ("{}".format(one_peak.new_points_pos[i]), end="\t", file=listfile)
                    print ("{}".format(one_peak.new_points_pos[0]), end="\t", file=listfile)
                else: 
                    for i in range(1,s_dim):
                        print ("{}".format(one_peak.peak_points_pos[i]), end="\t", file=listfile)
                    print ("{}".format(one_peak.peak_points_pos[0]), end="\t", file=listfile)
                if one_peak.peak_intens:
                    print (one_peak.peak_intens, end=" ", file=listfile) 
                if type_list=="new" and one_peak.is_center == False:
                    print ("to check", file=listfile)
                else: print (file=listfile)
    return



def current_time():
    curr_time = time.strftime("%H:%M:%S", time.localtime())
    return curr_time














"""                      MAIN PROGRAM                      """
if __name__ == "__main__":
    # print ("\n\nSTART")
    """File Reading"""
    
   
    print("\n=== Reading input files ===")
    # print ("File reading - start")
    spectrum = CSpectrum(args.filename,args.peak_list)
    # Spectra_dim, Nucl_name, Points_number, Tile_size, Spectrometer_fq, SW_size, SW_size_ppm, data_center, N_Tiles = read_ucsf_info(filename)
    # Peaks = read_peaklist(args.peak_list, Spectra_dim)                       # reading position of peak from peak list in Sparky format
    with open('{}/info.txt'.format(peak_list_dir_new), 'w') as txtfile:
        txtfile.write("File name: {}\n\n".format(args.filename))
    print("Peak list read: {} peaks".format(len(spectrum.peaks)))
    print("=== Reading input files finished ===")
    peaks = spectrum.peaks
    if args.OnlyPoints:
        peaks.calc_peak_points(spectrum, "Float")
        Print_Peak_List_points(peaks, spectrum.spectra_dimentionality, "orgin")
    else:
        peaks.calc_peak_points(spectrum)
        Print_Peak_List_points(peaks, spectrum.spectra_dimentionality, "orgin")

        """Noise calculation"""
        if UserPeakLevelFlag:
            spectrum.SetUpPeakLevel()            
            print_raport("Peak level starts from = {:.2e}".format(args.peak_level))
        else:
            print ("\n=== Noise calculation ===")
            spectrum.CalcNoise()
            print_raport("Average noise level = {:.2e}\nPeaks cutline = {:.2e}\n".format(spectrum.noise_level,spectrum.peak_level))
            print ("=== Noise calculation finished ===")


        if args.OnlyNoise:
            pass
        else:    
            """Center peaks and read intens"""
            print ("\n=== Peak centering and intensity reading ===")
            Peak_not_moved = 0
            Peak_moved = 0
            Peak_not_visible = 0
            Vector_set2 = [-1,0,1]
            peak_centering_info = []
            peak_height_list = []

            for indexpeak, one_peak in enumerate(Peaks):
                Try_position = one_peak.peak_points_pos
                Orgin_pos = one_peak.peak_points_pos
                try_pos=[]
                try_intens=[]
                while one_peak.is_center == False:
                    List_of_around_points, Points_intens, try_index = intens_aroun_peak(Try_position, Orgin_pos, filename, Spectra_dim, Tile_size, N_Tiles)
                    try_pos.append(deepcopy(Try_position))
                    try_intens.append(deepcopy(Points_intens[try_index]))
                    if Points_intens[try_index] > 0 :
                        maximum = max
                    else: 
                        maximum = min    
                    if abs(Points_intens[try_index]) > noise_level or args.noRemoveFlag==True:
                        if Points_intens.index(maximum(Points_intens)) != try_index:            # if current position of peak is not the highest: index of the highest point is not equal 0 (0 is index of current peak position)
                            Try_position = List_of_around_points[Points_intens.index(maximum(Points_intens))]
                            # print ("Peak:", indexpeak, "must move")
                            one_peak.was_moved = True
                        else:
                            one_peak.is_center = "maybe"
                            one_peak.peak_intens = Points_intens[try_index]
                            # print ("Peak:", indexpeak, "is in the highest place")
                            
                            if one_peak.was_moved == True:
                                one_peak.new_points_pos = Try_position
                                Peak_moved +=1
                            else: Peak_not_moved +=1
                    else: 
                        one_peak.is_visible = False
                        Peak_not_visible += 1 
                        # print ("Peak:", indexpeak, "is not visible")
                        break

                #check if peak is in the highest position
                List_of_around_points, Points_intens, try_index = intens_aroun_peak(Try_position, Orgin_pos, filename, Spectra_dim, Tile_size, N_Tiles)
                if Points_intens[try_index] > 0 :
                    maximum = max
                else: 
                    maximum = min    
                if Points_intens.index(maximum(Points_intens)) == try_index:
                    one_peak.is_center = True
                    peak_height_list.append(abs(Points_intens[try_index]))
                else: 
                    one_peak.is_center = False

                if abs(Points_intens[try_index]) < peak_level and args.noRemoveFlag==False:
                    one_peak.is_visible = False
                    Peak_not_visible += 1 

                # with open('{}/info.txt'.format(peak_list_dir_new), 'a') as txtfile:
                if one_peak.is_visible == False:
                    peak_centering_info.append("\n{}\t({})\t- is not visible".format(indexpeak, one_peak.descript))
                else:
                    peak_centering_info.append("\n{}\t({})".format(indexpeak, one_peak.descript))

                for i in range(len(try_pos)):
                    peak_centering_info.append("\t{}\t{}".format(try_pos[i], try_intens[i]))
                if one_peak.is_center == False:
                    peak_centering_info.append("                     is not in the highest position\t")

            average_peaks_level = round(sum(peak_height_list)/len(peak_height_list), 1)
            print_raport("Peak not moved = {}\nPeak moved = {}\nPeak not visible = {}".format(Peak_not_moved,Peak_moved,Peak_not_visible))
            print_raport("Average peaks height = {:.2e}".format(average_peaks_level))
            if UserPeakLevelFlag==False:
                signal_to_noise = average_peaks_level/noise_level
                print_raport("Signal to noise ratio = {:.2f}".format(signal_to_noise))
            print_raport("\n")
            for i in peak_centering_info:
                print_raport(i, in_terminal=False)

            print ("=== Peak centering and intensity reading finished ===")


            """Print peak list"""
            print ("\n=== Printing new peak lists ===")

            NewPeakList_ppm = Print_Peak_List_ppm(Peaks, Spectra_dim, Points_number, SW_size_ppm, data_center)
            NewPeakList_points = Print_Peak_List_points(Peaks, Spectra_dim, "new")
            # print ("\nPrint peak list - finish")
