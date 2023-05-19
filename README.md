# CCR_scripts   
Set of scripts for working with CCR spectra (2D-4D) including:
<!-- - [CCR\_scripts](#ccr_scripts) -->
- [CCR\_scripts](#ccr_scripts)
  - [read\_ucsf.py](#read_ucsfpy)
  - [calc\_CCR\_rate.py](#calc_ccr_ratepy)
  - [read\_point\_list\_and\_compere.py](#read_point_list_and_comperepy)

## read_ucsf.py    
read_ucsf script after reading ucsf file and peak list (in Sparky format), check peak intensity of peak and if it is not in the highest position - move it. After this script prints peak lists in ppm value and points value.  

<strong> For run script type in command line: </strong>    
```
    python3 read_ucsf.py [ucsf path] [peak list path]
```

**Additionaly options:**        
```
--np [num]              - to change number of points for calculate noise level: N^(spectra dimentionality + 1); normally is N = 10  
--plevel [num]          - if you know level when starting appear, add this with scientific numer notation e.g. 1e+7     
--norm (or --noRemove)    - add this if you do not want remove invisible peaks    
--onlypoints            - add this if you want only change ppm value to points value    
--o (or --output_name)    - add this if you want specific output name     
--noise                 - add this if you want calculate only noise level      
```

<br><br>

## calc_CCR_rate.py

The script to calculate CCR constants. In director has to contains:     
- peak lists (with peak position in points and peak height)     
- experiment_set        
- sequence in FASTA format      
<br>

<strong> For run script type in command line: </strong> 
```
python3 calc_CCR_rate.py [file director]
```

**Additionaly options:**    
```
    --seq       - if name of file with aminoacid sequence is not 'seq', add this with name of file       
    --refgamma  - if you have file with reference values of CCR rates, add this with name of file 
                (file must be .csv, columns name should be: AA, psi_angle (or/and phi_angle), CCR name )  
    --expset    - if you want use experiments setup file with diffrent filename than "experiments_set.txt" 
                (structure of file must be the same as orginal file)         
```

<br><br>

## read_point_list_and_compere.py

<strong> For run script type in command line: </strong>    
```
python3 read_point_list_and_compere [peak list1 path] [peak list2 path] [peak list3 path]...[peak listN path] --dim []
```
where:  
-dim [] - dimentionality of peak list

**Additionaly options:**    
```
    --comp2list [num] [num] - to additionaly compare two of peak list from bove       
    --out_dir [path]        - path to direction where output file will be put     
    --fl [path]             - path to file with list of relativ path to peaklists - every set of file must separate by "------"          
    --intens                - print peak height, no peak possition
```