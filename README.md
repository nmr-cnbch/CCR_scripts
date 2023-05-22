Table of contents:
- [CCR\_scripts](#ccr_scripts)
  - [read\_ucsf.py](#read_ucsfpy)
  - [calc\_CCR\_rate.py](#calc_ccr_ratepy)
    - [Build up of the experiment\_set file](#build-up-of-the-experiment_set-file)
    - [CCR names](#ccr-names)
  - [read\_point\_list\_and\_compere.py](#read_point_list_and_comperepy)

# CCR_scripts   
Set of python (version 3) scripts for working with NMR spectra (UCSF Sparky format), including:
  - read_ucsf.py
  - calc_CCR_rate.py
  - read_point_list_and_compere.py

software required:
  - python 3
  - numpy library
  - matplotlib library
  - scipy library


## read_ucsf.py    
read_ucsf script after reading ucsf file and peak list (in Sparky format), check peak intensity of peak and if it is not in the highest position - move it. After this script prints peak lists in ppm value and points value.  

<strong> For run script type in command line: </strong>    
```
    python3 read_ucsf.py [ucsf path] [peak list path]
```

**Additionaly options:**        
```
--np [num]                - to change number of points for calculate noise level: N^(spectra dimentionality + 1); normally is N = 10  
--plevel [num]            - if it is known a level when peaks are starting appear, add this with scientific numer notation e.g. 1e+7     
--norm (or --noRemove)    - without removing invisible peaks    
--onlypoints              - only change ppm value to points value    
--o (or --output_name) [file name]   - specific output name     
--noise                   - only calculate noise level      
```
File required:
  - multidimentional NMR spectra in ucsf-Sparky format
  - peak list in Sparky format  



<br><br>

## calc_CCR_rate.py

The script to calculate cross-correlated relaxation (CCR) rates, measure by quantidative approach (two spectra: reference and transfer). An eqation used for calculating CCR rates was:   
$\Gamma = \frac{1}{Tc}\cdot arctanh(\frac{I_{trans}}{I_{ref}})$ 

**File required:**  
- peak lists (with peak position in points and peak height)     
- experiment_set        
- sequence in FASTA format (seq file)     
  
<br>

<strong> For run script type in command line: </strong>   
```
    python3 calc_CCR_rate.py [file director]
```

**Additionaly options:**    
```
--seq [file name]      - if name of file with aminoacid sequence is not 'seq', add this with name of file       
--refgamma [file name] - compering calculated CCR rates with reference values  
          (file must be .csv, columns name should be: AA, psi_angle (or/and phi_angle), CCR name )  
--expset [file name]   - if you want use experiments setup file with diffrent filename than "experiments_set.txt" (structure of file must be the same as orginal file)         
```



### Build up of the experiment_set file
```
type_of_CCR: 'CCR_AA'         # name of CCR, this name will be use do make charts and peaklist
dir_auto: 10                  # dir name of auto version 
dir_cross: 25                 # dir name of cross version 
dimension: 4                  # dimensionality of spectra
nucl: N CO CA H               # nuclei of all peaks: H, N, C, CA, CB, HA, HB
pos_nucl: 1 -1 -2 1           # nuclei position of all peaks: -2, -1, 0, 1, 2 
angle_num:    2               # number of measured angle: 1, 2
angle_names: phi psi          # angle: phi, psi
angle_pos: 1 -1               # angle position of all peaks: -2, -1, 0, 1, 2
NS_auto: 3                    # number of scans in auto version
NS_cross: 10                  # number of scans in cross version
TC: 0.0286                    # time of ccr evolution 
H_roi: 6.0 9.6                # downfield and upfield of direct dimention
other: 100NUS                 # add this as a readable comment for the script to particular experiment
```

**Three types of comments**   
1. "other: [comment text]" - add this as a readable comment for the script to particular experiment.
2. "#" - inline comment, script will not read from this sign to the end of line   
3. "====" - block comment, script will not read from this sign to the end for file 


**Additionaly options:**  
If you want use specific name of peak list replace dir_auto and dir_cross with auto_name and cross_name and after ":" write name of peak list.  

example:
```
auto_name: 5_1
cross_name: 6_1
```
You can compere diffrent data for the same CCR rate. If the script finds 2 data sets for the same CCR rate then it will prepare a graph comparing them. Remember to add comments to specific experiments in the experiment_set file like:
``` 
other: [number - 3 digits][parameter - 3 letters]
``` 
example:  
```
other: 100NUS
```     

### CCR names   
Our scrips can work with any type of CCR rate but we builed up some of them:
- CCR_1 - H<sup>N</sup><sub>i</sub>N<sub>i</sub> DD – C<sup>$\alpha$</sup><sub>i-1</sub>H<sup>$\alpha$</sup><sub>i-1</sub> DD  
- CCR_2 - H<sup>N</sup><sub>i</sub>N<sub>i</sub> DD – C<sup>$\alpha$</sup><sub>i</sub>H<sup>$\alpha$</sup><sub>i</sub> DD  
- CCR_3 - H<sup>N</sup><sub>i</sub>N<sub>i</sub> DD - H<sup>N</sup><sub>i-1</sub>N<sub>i-1</sub> DD  
- CCR_4 - H<sup>N</sup><sub>i</sub>H<sup>$\alpha$</sup><sub>i-1</sub> DD – C'<sub>i-1</sub> CSA 
- CCR_5 - C<sup>$\alpha$</sup><sub>i-1</sub>H<sup>$\alpha$</sup><sub>i-1</sub> DD – C'<sub>i-1</sub> CSA  
- CCR_6 - C<sup>$\alpha$</sup><sub>i</sub>H<sup>$\alpha$</sup><sub>i</sub> DD – C'<sub>i-1</sub> CSA  
- CCR_7 - H<sup>N</sup><sub>i</sub>N<sub>i</sub> DD – C'<sub>i</sub> CSA  
- CCR_8 - C'<sub>i-1</sub> CSA –  C'<sub>i</sub> CSA   

Additionaly for CCR 5 and CCR 6 we multipy CCR rate by -1, which is consistent with the experiment.


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