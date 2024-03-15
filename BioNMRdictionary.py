import sys


def Nucl1to2(nucl:str) -> str:
    if nucl=="H": i="HN"
    elif nucl=="C": i="CO"
    
    elif nucl=="N": i="N"
    elif nucl=="CA": i="CA"
    elif nucl=="CB": i="CB"
    elif nucl=="HA": i="HA"
    elif nucl=="HB": i="HB"
    else:
        print ("Wrong name of nucleus (",nucl,")")
        sys.exit(1)
    return i


def Res1to3(res:str) -> str:
    if res=="P": i="PRO"
    elif res=="A": i="ALA"
    elif res=="S": i="SER"
    elif res=="T": i="THR"
    elif res=="C": i="CYS"
    elif res=="C": i="CYSox"
    elif res=="C": i="CYSred"
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
    
    elif res=="-": i=" "
    elif res=="": i=" "
    
    else:
        print ("Wrong name of residue (",res,"), check the seq file")
        sys.exit(1)
    return i

def Res3to1(res:str) -> str:
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

def CCRname2PrettyRateNamePLT(ccrname):
    if ccrname=="CCR_1": i="C$^{a}_{i-1}H^{a}_{i-1}$ - N$_{i}H^{N}_{i}$"          #
    elif ccrname=="CCR_2": i="C$^{a}_{i}H^{a}_{i}$ - N$_{i}H^{N}_{i}$" 
    elif ccrname=="CCR_3": i="N$_{i-1}H^{N}_{i-1}$ - N$_{i}H^{N}_{i}$"
    elif ccrname=="CCR_4": i="$H^{N}_{i}H^{a}_{i-1}$ - C$_{i-1}$"
    elif ccrname=="CCR_5": i="C$^{a}_{i}H^{a}_{i}$ - C$_{i}$"
    elif ccrname=="CCR_6": i="C$^{a}_{i}H^{a}_{i}$ - C$_{i-1}$"
    elif ccrname=="CCR_7": i="N$_{i}H^{N}_{i}$ - C$_{i}$"
    elif ccrname=="CCR_8": i="C$_{i-1}$ - C$_{i}$"
    elif ccrname=="CCR_9": i="C$^{a}_{i-1}H^{a}_{i-1}$ - C$^{a}_{i}H^{a}_{i}$"
    else:
        print ("Wrong name of CCR (",ccrname,")for matplotlib module")
        sys.exit(1)
    return i