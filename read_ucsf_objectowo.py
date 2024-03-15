#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
from turtle import circle

from more_itertools import first


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
import random

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
                    help="add this if you do not want remove invisible __peaks")

parser.add_argument("-op", "--onlypoints", dest='OnlyPoints', action="store_true",  default=False,
                    help="add this if you want only change ppm value to points value")

parser.add_argument("-o", "--output_name", type=Path, 
                    help="add this if you want specific output name")

parser.add_argument("-n", "--noise", dest='OnlyNoise', action="store_true",  default=False,
                    help="add this if you want calculate only noise level")

parser.add_argument("-sn", "--signal3noise", dest='SignalToNoise', type=float,
                    help="add this if you want setup minimal signal to noise ratio")

parser.add_argument("-rec", "--reconstrutedspecrum", dest='ReconstructionFlag', action="store_true",  default=False,
                    help="""add this if you your spectrum was reconstructed. 
                    For reconstracted spectra we calculate 'noise' by measure random points on H-cross-section with peaks. 
                    Otherwise, for traditional collect spectra noise is calculated by measure random points from across the spectrum""")

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
    print ("Peak centering without removing invisible __peaks")
if args.OnlyPoints:
    print ("OnlyPoints is on")    
if args.OnlyNoise:
    print ("Only calculate noise level option is on")    
if args.SignalToNoise:
    print (f"Minimal signal to noise ration setup to {args.SignalToNoise}") 

if args.ReconstructionFlag:
    noise_type="artifacts"
    print (f"Reconstructed spectrum - noise will calculate by measure random points on H-cross-section with peaks") 
else:
    noise_type="termal"



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
    __slots__ = [ '__filename', '__nucl_name','__points_num','__tile_size','__n_tiles',
                 '__spectr_fq','__sw','__sw_ppm','__data_center','__spectra_dim',
                 '__peaks','__peak_level','__noise_level','__minimal_signal_to_noise']

    def __init__(self,filename,peak_list,**kwargs):
        self.__filename = filename
        self.__nucl_name = []
        self.__points_num = []
        self.__spectra_dim = 0 #spectra_dimentionality
        self.__tile_size = []
        self.__n_tiles=[]      # number of tiles along axis [direct, 1, 2, 3] ,# point number with peak in tile  
        self.__spectr_fq =[]
        self.__sw = []
        self.__sw_ppm = []
        self.__data_center = []
        # self.__s_dim = 0  #spectra_dimentionality
        self.__peak_level = 0.0
        self.__noise_level = 0.0
        self.__minimal_signal_to_noise = 0.0
        
        self.read_specra_paramiters()
        self.calc_spectra_width_in_ppm()

        if "signal2noise" in kwargs:
            self.__minimal_signal_to_noise = kwargs["signal2noise"]
        
        self.__peaks = self.read_peaklist(peak_list)

         
    @property
    def spectra_dim(self) -> list:
        return self.__spectra_dim
    
    @property
    def tile_size(self) -> list:
        return self.__tile_size

    @property
    def n_tiles(self) -> int:
        return self.__n_tiles
        
    def read_specra_paramiters(self):
        """Method to read paramiters from UCSF file (binary file)

           Output: 
        """
        with open(self.__filename, "rb") as ucsf_file:
            """Read spectra dimentionality"""
            ucsf_file.seek(10)
            ucsf_data = ucsf_file.read(1)                                   # 1 byte per information 
            self.__spectra_dim = int.from_bytes(ucsf_data, byteorder='big')         # convert bytes to integer value 
            print ("Spectra dimensiolity:", self.__spectra_dim)

            """Read spectra format"""
            ucsf_file.seek(13)
            ucsf_data = ucsf_file.read(1)                                   # 1 byte per information 
            spectra_format = int.from_bytes(ucsf_data, byteorder='big')         # convert bytes to integer value  

            """Read nucleus name"""
            for i in range(0,self.__spectra_dim):
                ucsf_file.seek(180+128*i)                                       # nucleus name (1H, 13C, 15N, 31P); first is after 180 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(6)                                   # 6 byte per information 
                text_data = ucsf_data.decode('utf-8')                           # convert bytes to string    
                self.__nucl_name.append(deepcopy(text_data))
            self.__nucl_name.insert(0,self.__nucl_name.pop())
            print ("Nucleus names:",*self.__nucl_name)

            """Read number of points per asix"""
            for i in range(0,self.__spectra_dim):
                ucsf_file.seek(188+128*i)                                       # integer number of data points along axis; first is after 188 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(4)                                   # 4 byte per information 
                bytes2int = int.from_bytes(ucsf_data, byteorder='big')          # convert bytes to integer value  
                self.__points_num.append(deepcopy(bytes2int)) 
            self.__points_num.insert(0,self.__points_num.pop())
            print ("Number of points in axis:",*self.__points_num)

            """Read integer tile size along this axis"""
            for i in range(0,self.__spectra_dim):
                ucsf_file.seek(196+128*i)                                       # integer tile size along this axis; first is after 196 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(4)                                   # 4 byte per information 
                bytes2int = int.from_bytes(ucsf_data, byteorder='big')          # convert bytes to integer value  
                self.__tile_size.append(deepcopy(bytes2int))
            self.__tile_size.insert(0,self.__tile_size.pop())

            """Read float spectrometer frequency for this nucleus (MHz) """
            for i in range(0,self.__spectra_dim):
                ucsf_file.seek(200+128*i)                                       # float spectrometer frequency for this nucleus (MHz) ; first is after 196 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(4)                                   # 4 byte per information 
                [bytes2float] = struct.unpack('>f', ucsf_data)                  # convert bytes to float value
                self.__spectr_fq.append(deepcopy(bytes2float))
            self.__spectr_fq.insert(0,self.__spectr_fq.pop())
            print ("Spectrometer frequency (MHz): ",' '.join("{:.2f}".format(x) for x in self.__spectr_fq))

            """Read float spectral width (Hz)"""
            for i in range(0,self.__spectra_dim):
                ucsf_file.seek(204+128*i)                                       # float spectral width (Hz) ; first is after 196 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(4)                                   # 4 byte per information 
                [bytes2float] = struct.unpack('>f', ucsf_data)                  # convert bytes to float value
                self.__sw.append(deepcopy(bytes2float))
            self.__sw.insert(0,self.__sw.pop())
            print ("Spectral width (Hz):", ' '.join("{:.0f}".format(x) for x in self.__sw))


            """Read float center of data (ppm)"""
            for i in range(0,self.__spectra_dim):
                ucsf_file.seek(208+128*i)                                       # float center of data (ppm) ; first is after 196 bytes, next one is after additional 128 bytes
                ucsf_data = ucsf_file.read(4)                                   # 4 byte per information 
                [bytes2float] = struct.unpack('>f', ucsf_data)                  # convert bytes to float value
                self.__data_center.append(deepcopy(bytes2float))
            self.__data_center.insert(0,self.__data_center.pop())
   

    def calc_spectra_width_in_ppm(self) -> None:
        """Calculate spectral width from Hz to ppm

           Output: self.__sw_ppm, self.__n_tiles
        """
        for p in range(len(self.__sw)):
            sw_sparky = self.__sw[p]-(self.__sw[p]/self.__points_num[p])
            self.__sw_ppm.append(deepcopy(sw_sparky/self.__spectr_fq[p]))
        print ("Spectral width (ppm):", ' '.join("{:.1f}".format(x) for x in self.__sw_ppm))
        
        for i in range (0,self.__spectra_dim):
            if self.__points_num[i]%self.__tile_size[i]==0:
                self.__n_tiles.append(deepcopy(int(self.__points_num[i]/self.__tile_size[i])))
            else: 
                self.__n_tiles.append(deepcopy(int(self.__points_num[i]/self.__tile_size[i])+1))

    def read_peaklist(self, peak_list) -> list[CPeak]:
        with open(peak_list, 'r') as pl:
            p_lines = pl.readlines()
            p_list = []
            for indexl, line in enumerate(p_lines):
                if indexl > 1 :
                    p_pos = CPeak(line,self.__spectra_dim)
                    p_list.append(deepcopy(p_pos))
                    # print (p_pos.peak_ppm_pos)     
        print("Peak list read: {} peaks".format(len(p_list)))  
        return p_list
        
    def ppm2points(self, ppm:float, which_dim:int) -> int: 
        sw_div_fnz = self.__sw_ppm[which_dim]/(self.__points_num[which_dim]-1)
        downfield = self.__data_center[which_dim]+self.__points_num[which_dim]/2*sw_div_fnz
        point_value = round((downfield-ppm)/(sw_div_fnz))
        if point_value > self.__points_num[which_dim]:
            point_value = self.__points_num[which_dim] - point_value
        if point_value < 0:
            print(f"przed {point_value}, self.__points_num[which_dim]: {self.__points_num[which_dim]}, sw_div_fnz: {sw_div_fnz}, downfield: {downfield}, ppm: {ppm}")
            print(f"""sw_div_fnz = self.__sw_ppm[which_dim]/(self.__points_num[which_dim]-1):
            {sw_div_fnz} = {self.__sw_ppm[which_dim]}/({self.__points_num[which_dim]}-1)""")
            point_value = self.__points_num[which_dim] + point_value
            print("po ",point_value)
        return point_value

    def points2ppm(self,point_value, which_dim:int) -> float: 
        sw_div_fnz = self.__sw_ppm[which_dim]/(self.__points_num[which_dim]-1)
        downfield = self.__data_center[which_dim]+self.__points_num[which_dim]/2*sw_div_fnz
        ppm_value = downfield-point_value*sw_div_fnz
        return ppm_value

    def calc_peaks_positions_in_points(self):   
        for onepeak in self.__peaks:
             onepeak.calc_peak_points(self)


    def find_points_aroun_peak(self,try_position:list[int], orgin_pos:list[int],distance=1) -> tuple[list[list[int]],int]:
        """Method for prepering list of points around of peak
           and list index where is the position of the place that is being checked for being the top of the peak 
            
           Output: list of points around of peak (list[[dim1,dim2...]...]), list index of checking place
        """
        vector_set2 = [x for x in range(-distance,distance+1)]
        list_of_around_points = []
        for k in range(self.__spectra_dim):
            for i in vector_set2:
                TestList = try_position[:]
                TestList[k] += i
                if abs(TestList[k]-orgin_pos[k])<=distance*2 and TestList[k] not in list_of_around_points:
                    list_of_around_points.append(deepcopy(TestList))
        try_index = list_of_around_points.index(try_position)
        return list_of_around_points, try_index
    
    def dispatch_generate_vec(self,distance):
        if self.__spectra_dim == 2:
            return generate_2Dvec_set(distance)
        elif self.__spectra_dim == 3:
            return generate_3Dvec_set(distance)
        elif self.__spectra_dim == 4:
            return generate_4Dvec_set(distance)
        elif self.__spectra_dim == 5:
            return generate_5Dvec_set(distance)

    def find_points_in_circle(self,orgin_pos:list[int],distance=5) -> tuple[list[list[int]],list,int]:
        """Method for prepering list of points around of peak
           and list index where is the position of the place that is being checked for being the top of the peak 
            
           Output: list of points around of peak (list[[dim1,dim2...]...]), list index of checking place
        """
        list_of_around_points = [orgin_pos]
        circle_encounter = []
        for one_circle in range(1,distance+1):
            circle_encounter.append(deepcopy(len(list_of_around_points)))
            vector_set = self.dispatch_generate_vec(one_circle)
            #print_raport(f"spectra dim: {self.__spectra_dim} circle number: {one_circle}, vector_set: {vector_set}")
            for one_vec in vector_set:   #type: ignore
                Test_pos = orgin_pos[:]
                if len(one_vec) == len(Test_pos):
                    for k in range(self.__spectra_dim):
                        if Test_pos[k]+one_vec[k]<= self.__points_num[k]:
                            Test_pos[k] += one_vec[k]
                    if calc_distance(Test_pos,orgin_pos)<=(one_circle*math.sqrt(2)) and Test_pos not in list_of_around_points:
                        list_of_around_points.append(deepcopy(Test_pos))
        try_index = list_of_around_points.index(orgin_pos)
        if try_index != 0:
            print(f"Starting point is not at first place = {try_index}")
        return list_of_around_points, circle_encounter, try_index
    

    def find_points_in_peak_line(self,try_position:list[int], orgin_pos:list[int],distance=1) -> tuple[list[list[int]],int]:
        """Method for prepering list of points lin the line of peak throu all dimentions
           and list index where is the position of the place that is being checked for being the top of the peak 
            
           Output: list of points around of peak (list[[dim1,dim2...]...]), list index of checking place
        """
        # vector_set2 = [-1,0,1]
        vector_set2 = [x for x in range(-distance,distance+1)]
        list_of_around_points = []
        for k in range(self.__spectra_dim):
            for i in vector_set2:
                TestList = try_position[:]
                TestList[k] += i
                if abs(TestList[k]-orgin_pos[k])<=distance*2:
                    list_of_around_points.append(deepcopy(TestList))
        try_index = list_of_around_points.index(try_position)
        # points_intens = self.read_intens(list_of_around_points)
        return list_of_around_points, try_index


    def find_noise_point(self, number_of_points_for_noise) -> list[list[int]]:
        """Method to find N^(s_dim+1) random points in spectrum for calculation of noise averange level

           Output: list of points: list[[dim1,dim2...]...]
        """
        dumpster=[]                              # list of wrong points: e.g. it's peak or really close to peak
        for i in range(len(self.__peaks)):
            dumpster.append(deepcopy(self.__peaks[i].peak_points_pos))
        suspect_peaks=[]                                        # list of peaks which probably are too close to our random point
        list_of_points=[]                 # list of good points with noise

        first_try = [-1]*self.__spectra_dim

        print ("Number of points for noise:", number_of_points_for_noise**(1+self.__spectra_dim))
        ii=1
        while len(list_of_points) < number_of_points_for_noise**(1+self.__spectra_dim):
            for f in range(len(first_try)):
                first_try[f]=randint(0,self.__points_num[f]-1) # __points_num[f]-1 because it takes too big number
            # print (first_try)
            if first_try not in dumpster:
                # print ("first_try - ok")
                for indexi, one_peak in enumerate(self.__peaks):
                    if first_try[0] in range (one_peak.peak_points_pos[0]-10, one_peak.peak_points_pos[0]+10):
                        suspect_peaks.append(deepcopy(one_peak.peak_points_pos))
                        for i in range(1,len(one_peak.peak_points_pos)):
                            if first_try[i] not in range (one_peak.peak_points_pos[i]-5, one_peak.peak_points_pos[i]+5):
                                suspect_peaks.pop()
                                break
                # print ("suspect_peaks", ii, suspect_peaks)

                if len(suspect_peaks) == 0:
                    list_of_points.append(deepcopy(first_try))
                else:
                    dumpster.append(deepcopy(suspect_peaks))
                suspect_peaks.clear()
                ii+=1
        # print (list_of_points)
        # print ("ile wyrzuciles?", len(__peaks), "--->", len(dumpster))
        return list_of_points


    def find_termal_noise_point(self, number_of_points_for_noise) -> list[list[int]]:
        """Method to find N^(s_dim+1) random points in spectrum for calculation of termal noise averange level

        Output: list of points: list[[dim1,dim2...]...]
        """
        dumpster=[]                       # list of wrong points: e.g. it's peak or really close to peak
        for one_peak in self.__peaks:
            org_peak_point_pos = one_peak.peak_points_pos
            # dumpster.append(deepcopy(org_peak_point_pos))
            for f in range(0,self.__spectra_dim):
                near_peak_point_pos, index = self.find_points_aroun_peak(org_peak_point_pos,org_peak_point_pos,distance=5)
                for point_set in near_peak_point_pos:
                    dumpster.append(deepcopy(point_set))
                
        list_of_points=[]                 # list of good points with noise
        
        possible_H_points = [x for x in range(0,self.__points_num[0]+1)]
        for indexi, one_peak in enumerate(self.__peaks):
            for point_to_remove in range(one_peak.peak_points_pos[0]-5,one_peak.peak_points_pos[0]+6):
                if point_to_remove in possible_H_points:
                    possible_H_points.remove(point_to_remove)

        # print(len(dumpster),len(possible_H_points))

        for Hposs in possible_H_points:
            for i in dumpster:
                if i[0] == Hposs:
                    print(i)

        first_try = [-1]*self.__spectra_dim
        # first_try:list[int]
        great_number_of_points_for_noise = number_of_points_for_noise**(1+self.__spectra_dim)
        ii=1
        while len(list_of_points) < great_number_of_points_for_noise:
            for f in range(1,len(first_try)):
                first_try[f]=randint(0,self.__points_num[f]-1) # __points_num[f]-1 because it takes too big number
            first_try[0] = random.choice(possible_H_points)
            # print (first_try)
            if first_try not in dumpster:
                GoodPoingFlag = [True]*self.__spectra_dim
                for indexi, one_peak in enumerate(self.__peaks):
                    for i in range(1,len(one_peak.peak_points_pos)):
                        if first_try[i] in range (one_peak.peak_points_pos[i]-5, one_peak.peak_points_pos[i]+5):
                            GoodPoingFlag[i] = False
                if True not in GoodPoingFlag:
                    dumpster.append(deepcopy(first_try))
                    print(first_try, "fail")
                else:
                    list_of_points.append(deepcopy(first_try))
                    
                ii+=1
                print(f"Number of points for termal noise: {len(list_of_points)}/{great_number_of_points_for_noise}", end="\r")
        print(f"Number of points for termal noise: {len(list_of_points)}/{great_number_of_points_for_noise}")
        return list_of_points


    def find_artefact_noise_point(self, number_of_points_for_noise) -> list[list[int]]:
            """Method to find N^(s_dim+1) random points in spectrum for calculation of artefact noise averange level
               Chosen points are from 

               Output: list of points: list[[dim1,dim2...]...]
            """
            great_number_of_points_for_noise = number_of_points_for_noise**(1+self.__spectra_dim)

            dumpster=[]                       # list of wrong points: e.g. it's peak or really close to peak
            for one_peak in self.__peaks:
                org_peak_point_pos = one_peak.peak_points_pos
                dumpster.append(deepcopy(org_peak_point_pos))
                for f in range(0,self.__spectra_dim):
                    near_peak_point_pos = deepcopy(org_peak_point_pos)
                    for point_num in range(one_peak.peak_points_pos[f]-5,one_peak.peak_points_pos[f]+6):
                        near_peak_point_pos[f] = point_num
                        dumpster.append(deepcopy(near_peak_point_pos))
            list_of_points=[]                 # list of good points with noise

            possible_H_points_list = []
            possible_point_list = []
            for one_peak in self.__peaks:
                for possible_H_point in range(one_peak.peak_points_pos[0]-1,one_peak.peak_points_pos[0]+2):
                    possible_H_points_list.append(deepcopy(possible_H_point))

            first_try = [-1]*self.__spectra_dim
            # print(len(dumpster),len(possible_H_points_list))
            # print ("Number of points for artefact noise:", number_of_points_for_noise**(1+self.__spectra_dim))
            ii=1
            while len(list_of_points) < great_number_of_points_for_noise:
                for f in range(1,len(first_try)):
                    first_try[f]=randint(0,self.__points_num[f]-1) # __points_num[f]-1 because it takes too big number
                first_try[0] = random.choice(possible_H_points_list)
                if first_try not in dumpster:
                    GoodPoingFlag = [True]*self.__spectra_dim
                    for indexi, one_peak in enumerate(self.__peaks):
                        for i in range(1,len(one_peak.peak_points_pos)):
                            if first_try[i] in range (one_peak.peak_points_pos[i]-5, one_peak.peak_points_pos[i]+5):
                                GoodPoingFlag[i] = False
                    # print(first_try)
                    if True not in GoodPoingFlag:
                        dumpster.append(deepcopy(first_try))
                        print(first_try, "fail")
                    else:
                        list_of_points.append(deepcopy(first_try))
                    ii+=1
                    print(f"Number of points for artefact noise: {len(list_of_points)}/{great_number_of_points_for_noise}", end="\r")
            print(f"Number of points for artefact noise: {len(list_of_points)}/{great_number_of_points_for_noise}")
            return list_of_points


    def find_artefact_noise_point_around_one_peak(self, number_of_points_for_noise,peak_number) -> list[list[int]]:
            """Method to find 0.5N^(s_dim+1) random points around peak for calculation of artefact noise averange level for this particual peak
               Chosen points are from 

               Output: list of points: list[[dim1,dim2...]...]
            """
            great_number_of_points_for_noise = number_of_points_for_noise**(1+self.__spectra_dim)/2

            dumpster=[]                       # list of wrong points: e.g. it's peak or really close to peak
            this_peak = self.__peaks[peak_number]

            org_peak_point_pos = this_peak.peak_points_pos
            dumpster.append(deepcopy(org_peak_point_pos))
            for f in range(0,self.__spectra_dim):
                near_peak_point_pos = deepcopy(org_peak_point_pos)
                for point_num in range(this_peak.peak_points_pos[f]-5,this_peak.peak_points_pos[f]+6):
                    near_peak_point_pos[f] = point_num
                    dumpster.append(deepcopy(near_peak_point_pos))
            list_of_points=[]                 # list of good points with noise

            possible_H_points_list = []
            for possible_H_point in range(this_peak.peak_points_pos[0]-1,this_peak.peak_points_pos[0]+2):
                possible_H_points_list.append(deepcopy(possible_H_point))

            first_try = [-1]*self.__spectra_dim
            ii=1
            while len(list_of_points) < great_number_of_points_for_noise:
                for f in range(1,len(first_try)):
                    first_try[f]=randint(0,self.__points_num[f]-1) # __points_num[f]-1 because it takes too big number
                first_try[0] = random.choice(possible_H_points_list)
                if first_try not in dumpster:
                    GoodPoingFlag = [True]*self.__spectra_dim
                    for one_peak in self.__peaks:
                        for i in range(1,len(one_peak.peak_points_pos)):
                            if first_try[i] in range (one_peak.peak_points_pos[i]-5, one_peak.peak_points_pos[i]+5):
                                GoodPoingFlag[i] = False
                    # print(first_try)
                    if True not in GoodPoingFlag:
                        dumpster.append(deepcopy(first_try))
                        print(first_try, "fail")
                    else:
                        list_of_points.append(deepcopy(first_try))
                    ii+=1
                    print(f"Number of points for artefact noise: {len(list_of_points)}/{great_number_of_points_for_noise}", end="\r")
            print(f"Number of points for artefact noise: {len(list_of_points)}/{great_number_of_points_for_noise}")
            return list_of_points


    def calc_list_of_points(self,peak_point:list[list[int]]) -> tuple[list[list[int]],list[list[int]]]:      
        """Method for converting points in the spectrum to positions in tiles

           Output: L - list with numbers of tile, R - list with point number in each tile; 
                   every position in both lists correspond to one peak (position in spectrum)
        """
        L = []                                        # tile number with our peak -1
        R = []                                        # point number with peak in tile  
        for one_peak_point in peak_point:
            l = []
            r = []
            for p in range(len(one_peak_point)):
                l.append(deepcopy(one_peak_point[p]//self.__tile_size[p]))
                r.append(deepcopy(one_peak_point[p]%self.__tile_size[p]))
            R.append(deepcopy(r))
            L.append(deepcopy(l))
        return L,R


    def read_intens(self, peak_pos_list:list[list[int]]) -> list[int]:
        """Method for reading peak height from UCSF file
            
           Output: list of peaks height
        """
        L, R = self.calc_list_of_points(peak_pos_list)       # L_list: tile number with our peak -1; R_list: point number with peak in tile; Tiles_num: number of tiles along axis [direct, 1, 2, 3]
        intens_list = []
        peak_points_pos_list = []

        for k in range(len(peak_pos_list)):
            
            '''  Calculation starting point in ucsf file for peak "k"  '''
            ppos = [-1.0]*(self.__spectra_dim*2)
            for D in range(1,self.__spectra_dim):
                pp = 1
                rr = 1
                for d in range(D+1,self.__spectra_dim):
                    pp = pp * self.__n_tiles[d]                        # Calculating terms with quantities and sizes of tiles
                    rr = rr * self.__tile_size[d]                      # Calculating terms with position points in a given tiles
                pp = pp * L[k][D] * numpy.prod(self.__tile_size) * self.__n_tiles[0] #type:ignore
                rr = rr * R[k][D] * self.__tile_size[0]
                ppos[D-1] = pp
                ppos[D-1+self.__spectra_dim]=rr
            ppos[self.__spectra_dim-1] = L[k][0] * numpy.prod(self.__tile_size) #type:ignore
            ppos[self.__spectra_dim*2-1] = R[k][0]
            ppos_sum = 180+128*self.__spectra_dim+sum(ppos)*4
            peak_points_pos_list.append(deepcopy(ppos_sum))
        #     print(f"\t Number of points found in UCSF tilles: {k+1}/{len(peak_pos_list)}", end="\r")
        # print(f"\t Number of points found in UCSF tilles: {len(peak_points_pos_list)}/{len(peak_pos_list)}")
        '''  Reading value from ucsf file for peak "k"  '''
        with open(self.__filename, "rb") as ucsf_file:
            for indexk, k in enumerate(peak_points_pos_list):
                try:
                    ucsf_file.seek(k)                                       
                    ucsf_data = ucsf_file.read(4)
                    [bytes2float] = struct.unpack('>f', ucsf_data)
                    intens_list.append(deepcopy(bytes2float))
                except OSError:
                    print(f"Invalid argument: {k} - skiping this point: {peak_pos_list[indexk]}")
                    continue
                except struct.error:
                    print(f"Invalid argument: {k} - skiping this point: {peak_pos_list[indexk]}")
                    continue
                    # print("""Invalid argument: {} - skiping this point
                    # List of points deposit in: {}/List_of_reading_points""".format(k,peak_list_dir_new))
                    # list_path = f"{peak_list_dir_new}/List_of_reading_points.txt"
                    # if not os.path.exists(list_path):
                    #     for 
                # finally:
                #     continue
            #     print(f"\t Number of readings of point heights : {indexk+1}/{len(peak_points_pos_list)}", end="\r")
            # print(f"\t Number of readings of point heights: {len(peak_points_pos_list)}/{len(peak_pos_list)}")
        return intens_list


    def read_intens_point(self, peak_pos_list:list[PointInUCSF]):
        """Method for reading peak height from UCSF file
            
           Output: list of peaks height
        """
        with open(self.__filename, "rb") as ucsf_file:
            for indexk, k in enumerate(peak_pos_list):
                try:
                    ucsf_file.seek(k.point_in_file)                                       
                    ucsf_data = ucsf_file.read(4)
                    bytes2float = float(struct.unpack('>f', ucsf_data)[0])
                    # print(f"bytes2float: {bytes2float}, k.point_intensity: {k.point_intensity}")
                    k.point_intensity = deepcopy(bytes2float)
                except OSError:
                    print(f"Invalid argument: {k} - skiping this point: {peak_pos_list[indexk]}")
                    continue
                except struct.error:
                    print(f"Invalid argument: {k} - skiping this point: {peak_pos_list[indexk]}")
                    continue
 


    def set_up_peak_level(self):
        if UserPeakLevelFlag:
            print_raport("Peak level starts from = {:.2e}".format(args.peak_level))
            self.__peak_level = args.peak_level
            self.__noise_level = deepcopy(self.__peak_level/50)
     

    def CalcNoise(self,noise_type="",peak_number=-1) -> None:
        """Method for calculation averange noise levela and minimal peak level (50*noise level)
            
           Output: self.__noise_level, self.__peak_level or self.__peaks[peak_number].noise_around_peak
        """
        if noise_type == "termal":
            print ("\n=== Noise calculation ===")
            list_of_noise_points = self.find_termal_noise_point(args.Number_of_points_for_noise)
            if self.__minimal_signal_to_noise == 0.0: 
                self.__minimal_signal_to_noise = self.__spectra_dim * 5
        elif noise_type == "artifacts":
            print ("\n=== Noise calculation ===")
            list_of_noise_points = self.find_artefact_noise_point(args.Number_of_points_for_noise)
            if self.__minimal_signal_to_noise == 0.0: 
                self.__minimal_signal_to_noise = self.__spectra_dim * 2.5
        elif peak_number != -1:
            print("\n",peak_number)
            list_of_noise_points = self.find_artefact_noise_point_around_one_peak(args.Number_of_points_for_noise,peak_number)
        else:
            print ("\n=== Noise calculation ===")
            list_of_noise_points = self.find_noise_point(args.Number_of_points_for_noise)

        Points_intens = self.read_intens(list_of_noise_points)
        print ("Reading noise height is ready")
        average_noise_level = round(sum(Points_intens)/len(Points_intens), 1)
        sum_of_squares=0.0
        for val in Points_intens:
            sum_of_squares+=(val-average_noise_level)**2
        noise_val = math.sqrt(sum_of_squares/len(Points_intens))
        
        if peak_number != -1:
            self.__peaks[peak_number].noise_around_peak = noise_val
            print_raport(f"Average noise level for peak number {self.__peaks[peak_number].descript} is: {self.__peaks[peak_number].noise_around_peak:.2e}")
        else:
            self.__noise_level = noise_val
            self.__peak_level = self.__noise_level * self.__minimal_signal_to_noise
            noise_type_text = f" {noise_type} "
            if args.OnlyNoise:
                with open(f"{peak_list_dir_new}/noise.txt", 'w') as noise_file:
                    noise_file.write("Average{}noise level = {:.2e}\nPeaks cutline = {:.2e}\n".format(noise_type_text,self.__noise_level,self.__peak_level))
            print_raport("Average{}noise level = {:.2e}\nPeaks cutline = {:.2e}\n".format(noise_type_text,self.__noise_level,self.__peak_level))
            print ("=== Noise calculation finished ===")


    def try_centering_this_peak(self,one_peak:CPeak,Orgin_pos:list,checking_distance:int)->tuple[list[int],int,dict]:
        # checking_distance = 5
        circle_Flag = [False]*(checking_distance+1) #if in this circle is highter place
        checking_places = []

        List_of_around_points, circle_encounter, start_index = spectrum.find_points_in_circle(Orgin_pos, distance=checking_distance+1)
        print_raport(f"\n\n{one_peak.descript}\ncircle_encounter: {circle_encounter}")

        for point_pos in List_of_around_points:
            p = PointInUCSF(points_pos=point_pos,spectrum=self)
            checking_places.append(deepcopy(p))
        # print(f"\n\n{one_peak.descript}\nList_of_around_points: {List_of_around_points}\ncircle_encounter: {circle_encounter}\nstart_index: {start_index}")
        
        spectrum.read_intens_point(checking_places)
        #print_raport(f"starting point: {checking_places[0].point_pos}, {checking_places[0].point_intensity:.2e}\n")
        #print_raport(f"Len of checking_places list: {len(checking_places)}\n")

        if len(List_of_around_points) != len(checking_places):
            print(f"""Length of list of points and list of intens is diffrent!
                    List_of_around_points = {len(List_of_around_points)}
                    Points_intens = {len(checking_places)}""")

        suspect_set = {}

        circle_start = circle_encounter[0]
        circle_end = circle_encounter[1]
        print_raport(f"circle_start: {circle_start}, circle_end: {circle_end}")
        for indexn in range(circle_start,circle_end+1):
            if checking_places[indexn].point_intensity > checking_places[start_index].point_intensity:
                #print(f"indexn: {indexn}, pos: {checking_places[indexn].point_pos} intens: {checking_places[indexn].point_intensity:.2e}")
                suspect_set[indexn] = [checking_places[indexn].point_pos, checking_places[indexn].point_intensity]  #[position, height]
                circle_Flag[0] = True
        # for indexn, n_point in enumerate(checking_places[circle_start:circle_end+1]):
            # if n_point.point_intensity > checking_places[start_index].point_intensity:
            #     print(f"indexn: {indexn}, pos: {n_point.point_pos} intens: {n_point.point_intensity}")
            #     suspect_set[indexn] = [n_point.point_pos, n_point.point_intensity]  #[position, height]
            #     circle_Flag[circle_number] = True
        
        if len(suspect_set) > 1:
            for circle_number in range(1,checking_distance):
                circle_start = circle_encounter[circle_number]
                circle_end = circle_encounter[circle_number+1]
                for indexn in range(circle_start,circle_end+1):
                    if checking_places[indexn].point_intensity > checking_places[start_index].point_intensity:
                        suspekt_peaks = deepcopy([x for x in suspect_set])
                        for suspect_point in suspekt_peaks:
                            if calc_distance(suspect_set[suspect_point][0],checking_places[indexn].point_pos)<2:
                                if indexn not in suspect_set and checking_places[indexn].point_intensity > suspect_set[suspect_point][1]:
                                    suspect_set[indexn] = [checking_places[indexn].point_pos, checking_places[indexn].point_intensity]  #[position, height]
                                    circle_Flag[circle_number] = True
                # for indexn, n_point in enumerate(checking_places[circle_start:circle_end+1]):
                #     if n_point.point_intensity > checking_places[start_index].point_intensity:
                #         suspekt_peaks = deepcopy([x for x in suspect_set])
                #         for suspect_point in suspekt_peaks:
                #             if calc_distance(suspect_set[suspect_point][0],n_point.point_pos)<2:
                #                 if indexn not in suspect_set and n_point.point_intensity > suspect_set[suspect_point][1]:
                #                     print(f"indexn: {indexn}, pos: {n_point.point_pos} intens: {n_point.point_intensity}")
                #                     suspect_set[indexn] = [n_point.point_pos, n_point.point_intensity]  #[position, height]
                #                     circle_Flag[circle_number] = True

            distance_from_starting_point = 0
            farthest_point_intes = 0
            # print(f"suspect_set: {suspect_set}\ncircle_Flag: {circle_Flag}")
            print_raport(f"circle_Flag: {circle_Flag}")

            for suspect_point in suspect_set:
                current_dist = calc_distance(suspect_set[suspect_point][0],checking_places[start_index].point_pos)
                if abs(suspect_set[suspect_point][1])>abs(farthest_point_intes):   # current_dist > distance_from_starting_point and 
                    distance_from_starting_point = deepcopy(current_dist)
                    farthest_point_pos = deepcopy(suspect_set[suspect_point][0])
                    farthest_point_intes = deepcopy(suspect_set[suspect_point][1])

            if circle_Flag[-1] == False:
                one_peak.is_center = "yes"
                one_peak.peak_intens = farthest_point_intes
                if farthest_point_pos != checking_places[start_index].point_pos:
                    one_peak.was_moved = True
                    one_peak.new_points_pos = farthest_point_pos
                    one_peak.peak_intens = farthest_point_intes
                    print_raport(f"farthest_point_pos: {farthest_point_pos}, farthest_point_intes: {farthest_point_intes}")
                    return farthest_point_pos, farthest_point_intes, suspect_set
            
                else:
                    one_peak.was_moved = True
                    one_peak.is_center = "yes"
                    one_peak.peak_intens = checking_places[start_index].point_intensity
                    one_peak.new_points_pos = checking_places[start_index]
                    return checking_places[start_index], checking_places[start_index].point_intensity, suspect_set
            else:
                one_peak.is_center = "no"
            
        else:
            one_peak.is_center = "yes"
            one_peak.new_points_pos = checking_places[start_index]
            one_peak.peak_intens = checking_places[start_index].point_intensity
            return checking_places[start_index], checking_places[start_index].point_intensity, suspect_set


    def Center_peaks_new(self,noise_type="termal"):
        # peaks_moving_dict = {"Peak_not_moved":0,"Peak_moved"}
        Peak_not_moved = 0
        Peak_moved = 0
        Peak_not_visible = 0
        peak_centering_info = []
        peak_height_list = []
        noise_list = []

        for indexpeak, one_peak in enumerate(self.__peaks):
            Orgin_pos = one_peak.peak_points_pos
            New_pos, New_intens, path_dict = spectrum.try_centering_this_peak(one_peak,Orgin_pos, checking_distance=5)
            
            if one_peak.is_center == "yes":
                peak_height_list.append(New_intens)

            if one_peak.was_moved:
                Peak_moved +=1
            else: Peak_not_moved +=1
            
            if noise_type == "artifacts":
                self.CalcNoise(peak_number=indexpeak)
                noise_list.append(deepcopy(one_peak.noise_around_peak))
                if abs(New_intens) < one_peak.noise_around_peak and args.noRemoveFlag==False:
                    one_peak.is_visible = False
                    Peak_not_visible += 1 
            else:
                if abs(New_intens) < self.__peak_level and args.noRemoveFlag==False:
                    one_peak.is_visible = False
                    Peak_not_visible += 1 

            # preraring output to info.txt 
            if one_peak.is_visible == False:
                peak_centering_info.append("\n{}\t({})\t- is not visible".format(indexpeak, one_peak.descript))
            else:
                peak_centering_info.append("\n{}\t({})".format(indexpeak, one_peak.descript))
            
            if len(path_dict)>0:
                for suspect_point in path_dict:
                    peak_centering_info.append("\t{}\t{}".format(path_dict[suspect_point][0], path_dict[suspect_point][1]))
            
            peak_centering_info.append("\t{}\t{}".format(New_pos, New_intens))
            
            if one_peak.is_center == "no":
                peak_centering_info.append("                     is not in the highest position\t")

        average_peaks_level = round(sum(peak_height_list)/len(peak_height_list), 1)
        print_raport("Peak not moved = {}\nPeak moved = {}\nPeak not visible = {}".format(Peak_not_moved,Peak_moved,Peak_not_visible))
        print_raport("Average peaks height = {:.2e}".format(average_peaks_level))

        if noise_type == "artifacts":
            self.__noise_level = round(sum(noise_list)/len(noise_list), 1)
        if UserPeakLevelFlag==False:
            signal_to_noise = average_peaks_level/self.__noise_level
            print_raport(f"Minimal signal_to_noise: {self.__minimal_signal_to_noise}")
            print_raport("Signal to noise ratio = {:.2f}".format(signal_to_noise))
        print_raport("\n")
        for i in peak_centering_info:
            print_raport(i, in_terminal=False)


    def CalcNoise_around_peaks(self) -> None:
        """Method for calculation averange noise level around particular peak
            
           Output: self.__peaks[peak_number].noise_around_peak
        """
        noise_peak_list = f"{peak_list_dir_new}/{peak_list_name}_peaks_noise.list"
        for peak_index in range(len(self.__peaks)):
            self.CalcNoise(peak_number=peak_index)
        with open(noise_peak_list, 'w') as peak_nois_file:
            print (f"  Peak description   \tPeak uncertainty (noise value around peak)", file=peak_nois_file)
            for one_peak in self.__peaks:
                print (f"{one_peak.descript}\t{one_peak.noise_around_peak:.4e}", file=peak_nois_file)


    def Print_Peak_List_ppm(self):
        """Method for preparing peak list (txt file in Sparky format)  
            
           Output: peak list in ppm
        """
        new_peak_list = "{0}/{1}_new_ppm.list".format(peak_list_dir_new,peak_list_name)
        print ("new_peak_list", new_peak_list)
        max_lenth_discrip = 0
        for p in self.__peaks:
            if len(p.descript)>max_lenth_discrip:
                max_lenth_discrip=len(p.descript)
        with open(new_peak_list, 'w') as listfile:
            print ("\tAssignment", end="", file=listfile) 
            for i in range(self.__spectra_dim):
                print ("\tw{}".format(i+1), end="", file=listfile)
            print ("\tData Height\t\t\t\tCentering peak list prepare from{}\n".format(peak_list_name), file=listfile)
            for indexpeak, one_peak in enumerate(self.__peaks):
                if one_peak.is_visible == True:
                    print ("{:{sentence_len}}".format(one_peak.descript, sentence_len=max_lenth_discrip), end="\t", file=listfile)
                    if one_peak.was_moved == True and one_peak.is_center != "no":
                        for i in range(1,self.__spectra_dim):
                            ppm_position = self.points2ppm(one_peak.new_points_pos[i],i)
                            print ("{:.3f}".format(ppm_position), end="\t", file=listfile)
                        ppm_position = self.points2ppm(one_peak.new_points_pos[0],0)
                        print ("{:.3f}".format(ppm_position), end="\t", file=listfile)
                    else: 
                        for i in range(1,self.__spectra_dim):
                            print ("{:.3f}".format(one_peak.peak_ppm_pos[i]), end="\t", file=listfile)
                        print ("{:.3f}".format(one_peak.peak_ppm_pos[0]), end="\t", file=listfile)
                    print (one_peak.peak_intens, end=" ", file=listfile) 
                    
                    if one_peak.is_center == "no":
                        print ("to check", file=listfile)
                    else: print (file=listfile)
        return


    def Print_Peak_List_points(self, type_list):
        """Method for preparing peak list (txt file in Sparky format)  
            
           Output: peak list in points
        """
        new_peak_list = "{}/{}_{}_points.list".format(peak_list_dir_new,peak_list_name,type_list)

        
        with open(new_peak_list, 'w') as listfile:
            print ("\tAssignment", end="", file=listfile) 
            for i in range(self.__spectra_dim):
                print ("\tw{}".format(i+1), end="", file=listfile)
            if type_list == "orgin":
                print ("\tData Height\t\t\t\t Orgin peaklist in points\n", file=listfile)
                # print ("{} peak list".format(type_list), new_peak_list)
            else:
                print ("\tData Height\t\t\t\tCentering peak list prepare from{}".format(peak_list_name), file=listfile)
                print ("Average noise level = {:.2e}\n".format(self.__noise_level), file=listfile)
                print ("{} peak list".format(type_list), new_peak_list)
            for one_peak in self.__peaks:
                if one_peak.is_visible == True:
                    print ("{}".format(one_peak.descript), end="\t", file=listfile)
                    if one_peak.was_moved == True and one_peak.is_center != "no":
                        for i in range(1,self.__spectra_dim):
                            print ("{}".format(one_peak.new_points_pos[i]), end="\t", file=listfile)
                        print ("{}".format(one_peak.new_points_pos[0]), end="\t", file=listfile)
                    else: 
                        for i in range(1,self.__spectra_dim):
                            print ("{}".format(one_peak.peak_points_pos[i]), end="\t", file=listfile)
                        print ("{}".format(one_peak.peak_points_pos[0]), end="\t", file=listfile)
                    if one_peak.peak_intens:
                        print (one_peak.peak_intens, end=" ", file=listfile) 
                    if type_list=="new" and one_peak.is_center == "no":
                        print ("to check", file=listfile)
                    else: print (file=listfile)
        return






class CPeak:
    def __init__(self,input_line,s_dim):
        self.peak_ppm_pos = []            # chemical shifts for all nuclei of peak, length depends on dimentionality
        self.peak_points_pos = []     # peak position in points
        self.is_visible = True         # if peak height is 50 time biger than noise level = True
        self.peak_intens = 0            # peak height  
        self.descript = ""              
        self.is_center = "no"          # if peak is in the highest point = "yes", 
        self.was_moved = False          # if peak was moved = True
        self.new_points_pos = []         # if peak was moved there will be new position in points
        self.new_ppm_pos = []           # if peak was moved there will be new position in ppm
        self.noise_around_peak = 0.0 

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

    def calc_peak_points(self, spectrum:CSpectrum):
        for p in range(len(self.peak_ppm_pos)):
            one_dim = spectrum.ppm2points(self.peak_ppm_pos[p],p)
            self.peak_points_pos.append(deepcopy(one_dim))


class PointInUCSF(CSpectrum):
    def __init__(self,points_pos:list,spectrum:CSpectrum):
        self._points_pos = points_pos #type: list
        self._tile_num = []
        self._point_in_tile = []
        self._point_in_file = 0
        self.point_intensity = 0.0 

        # print(f"CSpectrum.tile_size: {spectrum.tile_size}\nCSpectrum.n_tiles: {spectrum.n_tiles}\nCSpectrum.spectra_dim: {spectrum.spectra_dim}")
        for p in range(len(self._points_pos)):
            self._tile_num.append(deepcopy(self._points_pos[p]//spectrum.tile_size[p]))
            self._point_in_tile.append(deepcopy(self._points_pos[p]%spectrum.tile_size[p]))

        '''  Calculation starting point in ucsf file for peak "k"  '''
        ppos = [-1.0]*(spectrum.spectra_dim*2)
        for D in range(1,spectrum.spectra_dim):
            pp = 1
            rr = 1
            for d in range(D+1,spectrum.spectra_dim):
                pp = pp * spectrum.n_tiles[d]                        # Calculating terms with quantities and sizes of tiles
                rr = rr * spectrum.tile_size[d]                      # Calculating terms with position points in a given tiles
            pp = pp * self._tile_num[D] * numpy.prod(spectrum.tile_size) * spectrum.n_tiles[0] #type:ignore
            rr = rr * self._point_in_tile[D] * spectrum.tile_size[0]
            ppos[D-1] = pp
            ppos[D-1+spectrum.spectra_dim] = rr
        ppos[spectrum.spectra_dim-1] = self._tile_num[0] * numpy.prod(spectrum.tile_size) #type:ignore
        ppos[spectrum.spectra_dim*2-1] = self._point_in_tile[0]
        self._point_in_file = 180 + 128*spectrum.spectra_dim + sum(ppos)*4

    @property
    def point_pos(self):
        return self._points_pos
    
    @property
    def point_in_file(self):
        return self._point_in_file

"""Other functions"""



def current_time():
    curr_time = time.strftime("%H:%M:%S", time.localtime())
    return curr_time


def generate_2Dvec_set(distance:int) -> list:
    vec_set = []
    for d_1 in range(-distance, distance+1):
        for d_2 in range(-distance, distance+1):
            vec = [d_1,d_2]
            vec_set.append(deepcopy(vec))
    return vec_set


def generate_3Dvec_set(distance:int) -> list:
    vec_set = []
    for d_1 in range(-distance, distance+1):
        for d_2 in range(-distance, distance+1):
            for d_3 in range(-distance, distance+1):
                vec = [d_1,d_2,d_3]
                vec_set.append(deepcopy(vec))
    return vec_set


def generate_4Dvec_set(distance:int) -> list:
    vec_set = []
    for d_1 in range(-distance, distance+1):
        for d_2 in range(-distance, distance+1):
            for d_3 in range(-distance, distance+1):
                for d_4 in range(-distance, distance+1):
                    vec = [d_1,d_2,d_3,d_4]
                    vec_set.append(deepcopy(vec))
    return vec_set


def generate_5Dvec_set(distance:int) -> list:
    vec_set = []
    for d_1 in range(-distance, distance+1):
        for d_2 in range(-distance, distance+1):
            for d_3 in range(-distance, distance+1):
                for d_4 in range(-distance, distance+1):
                    for d_5 in range(-distance, distance+1):
                        vec = [d_1,d_2,d_3,d_4,d_5]
                        vec_set.append(deepcopy(vec))
    return vec_set


def calc_distance(peak1,peak2):
        sum_of_squares = 0
        for i in range(len(peak1)):
            # print (peak1, "and", peak2)
            sum_of_squares += (peak1[i]-peak2[i])**2
        dis =  sum_of_squares**(1/len(peak1))
        return dis


"""                      MAIN PROGRAM                      """
if __name__ == "__main__":
    # print ("\n\nSTART")
    """File Reading"""
    
   
    print("\n=== Reading input files ===")
    # print ("File reading - start")

    if args.SignalToNoise:
        spectrum = CSpectrum(args.filename,args.peak_list,signal2noise=args.SignalToNoise)  
    else:
        spectrum = CSpectrum(args.filename,args.peak_list)
    # Spectra_dim, Nucl_name, Points_number, Tile_size, Spectrometer_fq, SW_size, SW_size_ppm, __data_center, N_Tiles = read_ucsf_info(filename)
    # Peaks = read_peaklist(args.peak_list, Spectra_dim)                       # reading position of peak from peak list in Sparky format
    with open('{}/info.txt'.format(peak_list_dir_new), 'w') as txtfile:
        txtfile.write("File name: {}\n\n".format(args.filename))
    print("=== Reading input files finished ===")
    spectrum.calc_peaks_positions_in_points()
    spectrum.Print_Peak_List_points("orgin")
    if args.OnlyPoints == False:

        """Noise calculation"""
        if UserPeakLevelFlag:
            spectrum.set_up_peak_level()            
            print_raport("Peak level starts from = {:.2e}".format(args.peak_level))
        else:
            if args.ReconstructionFlag and args.OnlyNoise:
                spectrum.CalcNoise(noise_type=noise_type)
            elif args.ReconstructionFlag == False:
                spectrum.CalcNoise(noise_type=noise_type)

        if args.OnlyNoise == False:

            """Center __peaks and read intens"""
            print ("\n=== Peak centering and intensity reading ===")
            # spectrum.Center_peaks()

            spectrum.Center_peaks_new(noise_type=noise_type)
            print ("=== Peak centering and intensity reading finished ===")


            """Print peak list"""
            print ("\n=== Printing new peak lists ===")

            NewPeakList_ppm = spectrum.Print_Peak_List_ppm()
            NewPeakList_points = spectrum.Print_Peak_List_points("new")
            print ("\n=== Printing new peak lists finished ===")
            
            # print ("\n=== Noise around peaks calculation ===")
            # spectrum.CalcNoise_around_peaks()
            # print ("\n=== Noise around peaks calculation finished ===")


        # TO DO
        # dodac sprawdzanie czy lista pików pasuje do widma
        # dodac przesuwanie listy wzgledem widma
        # dodac opcje wyboru przycinania pikow lub zawijania
