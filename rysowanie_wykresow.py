# Importing libraries
import matplotlib.pyplot as plt
from copy import deepcopy
import sys
import numpy as np
from scipy import stats


Dir = sys.argv[1]


class CExp_data():
    def __init__(self):
        self.name="None"
        self.angle=[]
        self.gamma_calc=[]
        self.gamma_exp=[]
        self.min_max = [+100,-100]


def read_exp_file(file_name,dir):
    with open("{}/{}".format(dir,file_name), 'r') as gl:
        g_lines = gl.readlines()
        g_list = CExp_data()
        for indexl, line in enumerate(g_lines):
            if indexl == 0 :
                g_list.name = line
            elif indexl > 1:
                item = line.split()
                if len(item)==3:
                    g_list.angle.append(deepcopy(float(item[0])))
                    g_list.gamma_calc.append(deepcopy(float(item[1])))
                    g_list.gamma_exp.append(deepcopy(float(item[2])))
                    if float(item[1]) < g_list.min_max[0]:
                        g_list.min_max[0]=float(item[1])
                    if float(item[1]) > g_list.min_max[1]:
                        g_list.min_max[1]=float(item[1])
                    if float(item[2]) < g_list.min_max[0]:
                        g_list.min_max[0]=float(item[2])
                    if float(item[2]) > g_list.min_max[1]:
                        g_list.min_max[1]=float(item[2])
                    
    return g_list

def read_exp_file2(file_name,dir):
    with open("{}/{}".format(dir,file_name), 'r') as gl:
        g_lines = gl.readlines()
        g_list = CExp_data()
        for indexl, line in enumerate(g_lines):
            if indexl == 0 :
                g_list.name = line
            elif indexl > 1:
                item = line.split()
                if len(item)==2:
                    g_list.gamma_calc.append(deepcopy(float(item[0])))
                    g_list.gamma_exp.append(deepcopy(float(item[1])))
                    if float(item[0]) < g_list.min_max[0]:
                        g_list.min_max[0]=float(item[0])
                    if float(item[0]) > g_list.min_max[1]:
                        g_list.min_max[1]=float(item[0])
                    if float(item[1]) < g_list.min_max[0]:
                        g_list.min_max[0]=float(item[1])
                    if float(item[1]) > g_list.min_max[1]:
                        g_list.min_max[1]=float(item[1])
    return g_list

def read_theor_file(file_name,dir):
    with open("{}/{}".format(dir,file_name), 'r') as gl:
        g_lines = gl.readlines()
        g_list = CExp_data()
        for indexl, line in enumerate(g_lines):
            if indexl == 0 :
                g_list.name = line
            elif indexl > 1:
                item = line.split()
                if len(item)>0:
                    g_list.angle.append(deepcopy(float(item[0])))
                    g_list.gamma_calc.append(deepcopy(float(item[1])))
        
    return g_list

def LRegression_expresion(x,y):
    slope, intercept, r, p, std_err = stats.linregress(x, y)
    r2 = r**2
    lr_exp = '{:.3f}*x + {:.3f}'.format(slope,intercept)
    lr_exp_y_list = []
    for i in x:
        lr_exp_y_list.append(deepcopy(slope*i+intercept))
    return lr_exp_y_list, lr_exp, r2


if __name__ == "__main__":
    
    file_names = ["CCR_1","CCR_2","CCR_4","CCR_5","CCR_6"]
    angle_names_in_order = ['psi','phi','psi','psi','phi']
    file_names2 = ["CCR_1","CCR_2","CCR_3","CCR_4","CCR_5","CCR_6","CCR_7"]
    file_names_theor = []
    for fn in file_names:
        file_names_theor.append(deepcopy(fn+"_theor"))
    
    Tab = []

    for indexf, fn in enumerate(file_names):
        
        tab_calc = read_exp_file(fn,Dir)
        tab_theor = read_theor_file(file_names_theor[indexf],Dir)
        Tab.append(deepcopy([tab_calc,tab_theor]))

        # Plotting both the curves simultaneously
        plt.scatter(tab_calc.angle, tab_calc.gamma_exp, s=10, color='gray', label='experimental')
        plt.plot(tab_theor.angle, tab_theor.gamma_calc, color='black', label='structure-predicted')
        
        # Naming the x-axis, y-axis and the whole graph https://pythonforundergradengineers.com/unicode-characters-in-python.html
        if angle_names_in_order[indexf]=='psi':
            plt.xlabel('\u03A8, deg')
        elif angle_names_in_order[indexf]=='phi':
            plt.xlabel('\u03A6, deg')
        plt.ylabel('\u0393, $s^{-1}$')
        
        # Adding legend, which helps us recognize the curve according to it's color
        plt.legend(fontsize="small")
        
        # To load the display window
        plt.savefig("{}.png".format(fn), bbox_inches="tight", pad_inches=0.3, transparent=True)
        # plt.show()
        plt.clf()
        if fn == 'CCR_6':
            print ("mamy to")
            for indexi, i in enumerate(tab_calc.gamma_calc):
                print (tab_calc.angle[indexi], i, tab_calc.gamma_exp[indexi])
    
    Tab.append(deepcopy([read_exp_file2("CCR_3",Dir)]))
    Tab.append(deepcopy([read_exp_file2("CCR_7",Dir)]))

    for indextf, tfn in enumerate(Tab):
        
        # Plotting both the curves simultaneously
        plt.scatter(tfn[0].gamma_calc, tfn[0].gamma_exp, s=10, color='gray', label='structure-predicted vs experimental')
        plt.axline([0,0],slope=1, linestyle=(0, (5, 5)), linewidth=1.5, color='black', label='x=y')
        
        # Calculating Linear Regression
        lr_exp_y_list,lr_expresion, coefition, = LRegression_expresion(tfn[0].gamma_calc,tfn[0].gamma_exp)
        plt.plot(tfn[0].gamma_calc, lr_exp_y_list, color='cornflowerblue', label='{}, $R^2$ = {:.3f}'.format(str(lr_expresion),coefition))

        

        # Naming the x-axis, y-axis and the whole graph 
        plt.xlabel('gamma structure-predicted')
        plt.ylabel('gamma experimental')
        plt.xlim(tfn[0].min_max[0]-5,tfn[0].min_max[1]+5)
        plt.ylim(tfn[0].min_max[0]-5,tfn[0].min_max[1]+5)
        # plt.axis('square')
        plt.gca().set_aspect('equal', adjustable='box')
        
        # Adding legend, which helps us recognize the curve according to it's color
        plt.legend(fontsize="small")
        
        

        # To load the display window
        plt.savefig("{}_exp_vs_calc.png".format(tfn[0].name), bbox_inches="tight", pad_inches=0.3, transparent=True)
        # plt.show()
        plt.clf()

        
