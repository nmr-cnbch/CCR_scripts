#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
from ../BioNMRdictionary.py import *
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


filename = args.filename
if not os.path.exists(filename): 
    print("There is not such file like this: ", filename)
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

class Spectrum:
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
        self.peaks = read_peaklist(peak_list)

    
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
                        p_pos = CPeak()
                        p_list.append(deepcopy(p_pos))
                        # print (p_pos.peak_ppm_pos)       
            return p_list


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

        if type(input_line) == string:
            input_line = input_line.split()
        if type(input_line) == list:
            self.descript = item[0]
            # print ("p_pos.peak_ppm_pos", p_pos.peak_ppm_pos)
            self.peak_ppm_pos.append(deepcopy(float(item[s_dim])))
            for i in range(1,s_dim):
                self.peak_ppm_pos.append(deepcopy(float(item[i])))
        else:
            pass
            # raise ERROR 

    def ppm2points(spectrum:Spectrum, ppm:float, output_format:string): 
        sw_div_fnz = spectrum.sw_ppm/(points_num-1)
        downfield = spectrum.data_center+spectrum.points_num/2*sw_div_fnz
        if output_format == "Float":
            point_value = round((downfield-ppm)/(sw_div_fnz),2)
        else:
        point_value = round((downfield-ppm)/(sw_div_fnz))
        if point_value > spectrum.points_num:
            point_value = spectrum.points_num - point_value
        if point_value < 0:
            point_value = spectrum.points_num + point_value
        
        return point_value


    def points2ppm(point_value, spectrum:Spectrum): 
        sw_div_fnz = spectrum.sw_ppm/(spectrum.points_num-1)
        downfield = spectrum.data_center+spectrum.points_num/2*sw_div_fnz
        ppm_value = downfield-point_value*sw_div_fnz
        return ppm_value


    def calc_peak_points(points_num, peaks, sw_ppm, data_center,s_dim,num_format="Integer"):
        for p in range(len(peaks[peak_num].peak_ppm_pos)):
            one_dim = ppm2points(peaks[peak_num].peak_ppm_pos[p],points_num[p],data_center[p],sw_ppm[p],num_format)
            peaks[peak_num].peak_points_pos.append(deepcopy(one_dim))
        Print_Peak_List_points(peaks, s_dim, "orgin")
        return peaks


"""Reading functions"""










"""Calculating functions"""





def find_noise_point(peaks, s_dim, points_num, number_of_points_for_noise):

    dumpster=[]                              # list of wrong points: e.g. it's peak or really close to peak
    for i in range(len(peaks)):
        dumpster.append(deepcopy(peaks[i].peak_points_pos))
    suspect_peaks=[]                                        # list of peaks which probably are too close to our random point
    list_of_points=[]                 # list of good points with noise

    first_try = [None]*s_dim


    print ("Number of points for noise:", number_of_points_for_noise**(1+s_dim))
    ii=1
    while len(list_of_points)< number_of_points_for_noise**(1+s_dim):
        
        for f in range(len(first_try)):
            first_try[f]=randint(0,points_num[f]-1) # points_num[f]-1 because it takes too big number
        
        # print (first_try)
        if first_try not in dumpster:
            # print ("first_try - ok")
            for indexi, one_peak in enumerate(peaks):
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
    # print ("ile wyrzuciles?", len(peaks), "--->", len(dumpster))

    return list_of_points


def calc_list_of_points(peak_point, tile_size):
    L=[]                                        # tile number with our peak -1
    R=[]                                        # point number with peak in tile  
    for one_peak_point in peak_point:
        l=[]
        r=[]
        for p in range(len(one_peak_point)):
            l.append(deepcopy(one_peak_point[p]//tile_size[p]))
            r.append(deepcopy(one_peak_point[p]%tile_size[p]))
        R.append(deepcopy(r))
        L.append(deepcopy(l))

    return L, R


def read_intens(peak_pos_list, filename, s_dim, tile_size, n_tiles):

    L, R = calc_list_of_points(peak_pos_list, tile_size)       # L_list: tile number with our peak -1; R_list: point number with peak in tile; Tiles_num: number of tiles along axis [direct, 1, 2, 3]
    intens_list=[]
    peak_points_pos_list = []

    for k in range(len(peak_pos_list)):
        
        '''  Calculation starting point in ucsf file for peak "k"  '''

        ppos = [None]*s_dim*2
        for D in range(1,s_dim):
            pp = 1
            rr = 1
            for d in range(D+1,s_dim):
                pp = pp * n_tiles[d]                        # Calculating terms with quantities and sizes of tiles
                rr = rr * tile_size[d]                      # Calculating terms with position points in a given tiles
            pp = pp * L[k][D] * numpy.prod(tile_size) * n_tiles[0]
            rr = rr * R[k][D] * tile_size[0]
            ppos[D-1]=pp
            ppos[D-1+s_dim]=rr
        ppos[s_dim-1]=L[k][0] * numpy.prod(tile_size)
        ppos[s_dim*2-1]=R[k][0]
        ppos_sum=180+128*s_dim+sum(ppos)*4
        peak_points_pos_list.append(deepcopy(ppos_sum))

    '''  Reading value from ucsf file for peak "k"  '''
    with open(filename, "rb") as ucsf_file:
        for k in peak_points_pos_list:
            try:
                ucsf_file.seek(k)                                       
                ucsf_data = ucsf_file.read(4)
                [bytes2float] = struct.unpack('>f', ucsf_data)
                intens_list.append(deepcopy(bytes2float))
            except OSError:
                print("Invalid argument: {} - check points lists".format(k))

        

    return intens_list






def find_points_around(try_position, orgin_pos, vector_set, s_dim):
    list_of_points=[]

    for i in range (len(vector_set)):
        Flag = False
        one_point = []
        for j in range (s_dim):
            new_point = try_position[j] + vector_set[i][j]
            if abs(new_point-orgin_pos[j])<=2:
                one_point.append(deepcopy(new_point))
            else:
                Flag = True
        if Flag==False:
            list_of_points.append(deepcopy(one_point))


    return list_of_points


def intens_aroun_peak(try_position, orgin_pos, filename, s_dim, tile_size, n_tiles):
    vector_set2 = [-1,0,1]
    list_of_around_points = []
    for k in range(s_dim):
        for i in vector_set2:
            TestList = try_position[:]
            TestList[k] += i
            if abs(TestList[k]-orgin_pos[k])<=2:
                list_of_around_points.append(deepcopy(TestList))
    try_index = list_of_around_points.index(try_position)
    points_intens = read_intens(list_of_around_points, filename, s_dim, tile_size, n_tiles)

    return list_of_around_points, points_intens, try_index



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
    # item = peak_list.split("/")
    # peak_list_name = item[len(item)-1][:-5]
    # peak_list_dir_new = peak_list[:-(len(peak_list_name)+5)]+peak_list_name
    # if not os.path.exists(peak_list_dir_new):
    #     os.mkdir(peak_list_dir_new)
    new_peak_list = "{0}/{1}_{2}_points.list".format(peak_list_dir_new,peak_list_name,type_list)

    
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


# def calc_noise(peaks: CPeak, spectra_dim, points_number, number_of_points_for_noise, filename, tile_size:list, n_tiles:list):
#     list_of_noise_points = find_noise_point(peaks, spectra_dim, points_number, number_of_points_for_noise)
#     Points_intens = read_intens(list_of_noise_points, filename, spectra_dim, tile_size, n_tiles)
#     print ("list of points of noise is ready: {} points".format(len(Points_intens)))
#     average_noise_level = round(sum(Points_intens)/len(Points_intens), 1)
#     sum_of_squares=0.0
#     for val in Points_intens:
#         sum_of_squares+=(val-average_noise_level)**2
#     standard_deviation = math.sqrt(sum_of_squares/len(Points_intens))
#     return standard_deviation

# def calc_noise(peaks=Peaks, spectra_dim=Spectra_dim, points_number=Points_number, number_of_points_for_noise=args.Number_of_points_for_noise, 
#                filename=filename, tile_size=Tile_size, n_tiles=N_Tiles):

#     list_of_noise_points = find_noise_point(peaks, spectra_dim, points_number, number_of_points_for_noise)
#     Points_intens = read_intens(list_of_noise_points, filename, spectra_dim, tile_size, n_tiles)
#     # list_of_noise_points = find_noise_point(Peaks, Spectra_dim, Points_number, args.Number_of_points_for_noise)
#     # Points_intens = read_intens(list_of_noise_points, filename, Spectra_dim, Tile_size, N_Tiles)
#     print ("list of points of noise is ready: {} points".format(len(Points_intens)))
#     average_noise_level = round(sum(Points_intens)/len(Points_intens), 1)
#     sum_of_squares=0.0
#     for val in Points_intens:
#         sum_of_squares+=(val-average_noise_level)**2
#     standard_deviation = math.sqrt(sum_of_squares/len(Points_intens))
#     return standard_deviation


def calc_noise():
    list_of_noise_points = find_noise_point(Peaks, Spectra_dim, Points_number, args.Number_of_points_for_noise)
    Points_intens = read_intens(list_of_noise_points, filename, Spectra_dim, Tile_size, N_Tiles)
    print ("list of points of noise is ready: {} points".format(len(Points_intens)))
    average_noise_level = round(sum(Points_intens)/len(Points_intens), 1)
    sum_of_squares=0.0
    for val in Points_intens:
        sum_of_squares+=(val-average_noise_level)**2
    standard_deviation = math.sqrt(sum_of_squares/len(Points_intens))
    return standard_deviation









"""                      MAIN PROGRAM                      """
if __name__ == "__main__":
    # print ("\n\nSTART")
    """File Reading"""
    
   
    print("\n=== Reading input files ===")
    # print ("File reading - start")
    Spectra_dim, Nucl_name, Points_number, Tile_size, Spectrometer_fq, SW_size, SW_size_ppm, data_center, N_Tiles = read_ucsf_info(filename)
    Peaks = read_peaklist(args.peak_list, Spectra_dim)                       # reading position of peak from peak list in Sparky format
    with open('{}/info.txt'.format(peak_list_dir_new), 'w') as txtfile:
        txtfile.write("File name: {}\n\n".format(filename))
    print("Peak list read: {} peaks".format(len(Peaks)))
    print("=== Reading input files finished ===")

    if args.OnlyPoints:
        Peaks = calc_peak_points(Points_number, Peaks, SW_size_ppm, data_center, Spectra_dim, "Float")
        # OldPeakList_points = Print_Peak_List_points(Peaks, Spectra_dim, "orgin")
    else:
        Peaks = calc_peak_points(Points_number, Peaks, SW_size_ppm, data_center, Spectra_dim)
        OldPeakList_points = Print_Peak_List_points(Peaks, Spectra_dim, "orgin")

        """Noise calculation"""
        if UserPeakLevelFlag:
            print_raport("Peak level starts from = {:.2e}".format(args.peak_level))
            peak_level = args.peak_level
            noise_level = peak_level/50
        else:
            print ("\n=== Noise calculation ===")
            noise_level = calc_noise()
            peak_level = 50*noise_level
            print_raport("Average noise level = {:.2e}\nPeaks cutline = {:.2e}\n".format(noise_level,peak_level))
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

        #### Check intensity of points around the peak
        # TestList = [30, 31, 157, 107]
        # moves = [-2,-1,0,1,2]
        # with open('{}/intensmap.txt'.format(peak_list_dir_new), 'w') as txtmap:
        #     print ("Maps of intensity of points around peak:{}\n\n".format(TestList),file=txtmap)

        #     duzalista = []
        #     for k in range(4):
                
        #         malalista = []
        #         for i in moves:
        #             TestList = [30, 31, 157, 107]
        #             TestList[k] += i
        #             malalista.append(deepcopy(TestList))
        #         jedna_linia = read_intens(malalista, filename, Spectra_dim, Tile_size, N_Tiles)
        #         duzalista.append(deepcopy(malalista))
        #         print ("index = {}".format(k), end="\t", file=txtmap)
        #         for p in jedna_linia:
        #             print ("{:^20}".format(p), end="\t", file=txtmap)
        #         print ("\n\t\t", end="\t", file=txtmap)
        #         for b in malalista:
        #             print ("{:^20}".format(str(b)), end="\t", file=txtmap)
        #         print ("\n",file=txtmap)


        
        # print ("\n\nKONIEC\n\n")
        # Try_point_list=[692]
        # try_point(Try_point_list)