# CCR_scripts   
Set of Python scripts for working with NMR spectra (UCSF [Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format), including:
  - read_ucsf.py
  - calc_CCR_rate.py


### software required:
  - [Python 3.6](https://www.python.org/) with libraries:
    - [NumPy](https://www.numpy.org/) 
    - [Matplotlib](https://matplotlib.org/)
    - [SciPy](https://www.scipy.org/)



--- 
---  


## read_ucsf.py    
The script after read ucsf file and peak list (in [Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format), then check peak hight and if it is not in the highest position - move it. After this script prints peak lists in ppm value and points value.

### Launch in command line:  
```bash
python3 read_ucsf [-h] [-np NUMBER_OF_POINTS_FOR_NOISE] [-pl PEAK_LEVEL] [-nrm] [-op] [-o OUTPUT_NAME] [-n] [-sn SIGNALTONOISE] [-rec] ucsf_path peak_list_path
```

### Positional arguments:        
```bash
  ucsf_path             path to UCSF file
  peak_list_path        path to peak list in SPARKY format
```

### Optional arguments:      
```bash
  -h, --help            show this help message and exit
  -np, --npoints        to change number of points for calculate noise level: N^(spectra dimentionality + 1); normally is N = 10
  -pl, --peaklevel      if you know level when starting appear, add this with scientific numer notation e.g. 1e+7
  -nrm, --noRemove      add this if you do NOT want remove invisible peaks
  -op, --onlypoints     add this if you want only change ppm value to points value
  -o, --output_name     add this if you want specific output name
  -n, --noise           add this if you want calculate only noise level
  -sn, --signal3noise   add this if you want setup minimal signal to noise ratio
  -rec, --reconstrutedspecrum
                        add this if you your spectrum was reconstructed. For reconstracted spectra we calculate 'noise' by measure random points on H-cross-section with peaks. Otherwise, for traditional collect spectra noise is calculated by measure random points from across the spectrum
```
### File required:
  - multidimensional NMR spectra in ucsf-Sparky format
  - peak list in Sparky format  

### Output:
The directory with the following files:
- orginal peak list with converted chemical shifts to points in the spectrum (`name+"_orgin_points.list"`)
- two peak list after moving peaks, one contains peaks positions in chemical shifts (`name+"_new_ppm.list"`), another one with points in the spectrum (`name+"_new_points.list"`)
- peak list contains only peak names and uncertainity (`name+"_peaks_noise.list"`)
- text file (`info.txt`) with whole terminal output and additional informations to evaluate script work



<br><br>

## calc_CCR_rate.py

The script to calculate cross-correlated relaxation (CCR) rates, measure by quantitative approach (two spectra: reference and transfer). 
   
Our scrips can work with any type of CCR rate but we builed up some of them:
|   |   |
|---|---|     
| CCR_1 | H<sup>N</sup><sub>i</sub>N<sub>i</sub> DD – C<sup>α</sup><sub>i-1</sub>H<sup>α</sup><sub>i-1</sub> DD |      
| CCR_2 | H<sup>N</sup><sub>i</sub>N<sub>i</sub> DD – C<sup>α</sup><sub>i</sub>H<sup>α</sup><sub>i</sub> DD |      
| CCR_3 | H<sup>N</sup><sub>i</sub>N<sub>i</sub> DD - H<sup>N</sup><sub>i-1</sub>N<sub>i-1</sub> DD |     
| CCR_4 | H<sup>N</sup><sub>i</sub>H<sup>α</sup><sub>i-1</sub> DD – C'<sub>i-1</sub> CSA |    
| CCR_5 | C<sup>α</sup><sub>i-1</sub>H<sup>α</sup><sub>i-1</sub> DD – C'<sub>i-1</sub> CSA |    
| CCR_6 | C<sup>α</sup><sub>i</sub>H<sup>α</sup><sub>i</sub> DD – C'<sub>i-1</sub> CSA |    
| CCR_7 | H<sup>N</sup><sub>i</sub>N<sub>i</sub> DD – C'<sub>i</sub> CSA |     
| CCR_8 | C'<sub>i-1</sub> CSA –  C'<sub>i</sub> CSA |     


An equation used for calculating CCR rates was:   
$\Gamma = \frac{1}{Tc}\cdot arctanh(\frac{I_{trans}}{I_{ref}})$    
*Additionally for CCR 5 and CCR 6 we multiply CCR rate by -1, which is consistent with the experiment.*


### File required:  
- peak lists with peak position in points (`name+"_points.list"`) or chemical shifts (`name+"_ppm.list"`) and peak height 
- JSON file with experiments description   
- sequence in FASTA format (`seq` file)     

### Additional File: 
- peak lists with peak names and uncertainity (`name+"_peaks_noise.list"`)
- file with reference values of CCR rates

<br>

### Launch in command line:    
```bash
python3 calc_CCR_rate [-h] [-s SEQ_FILE_NAME] [-r REFGAMMA] [-e EXPSET] [-pub] [-pres] file_director
```

### Positional arguments:       
```bash
file_director         path to to director with all required
```

### Optional arguments:       
```bash
-h, --help              show this help message and exit
-s, --seq               if name of file with amino acid sequence is not 'seq', add this with name of file
-r, --refgamma          if you have file with reference values of CCR rates, add this with name of file (file must be .csv, columns name should be: AA, psi_angle (or/and phi_angle), CCR name)
-e, --expset            if you want use experiments setup file with different filename than /input.json/ (structure of file must be the same as orginal file)
-pub, --publication     if you want output in publication size
-pres, --presentation   if you want output in presentation size
```



### Build up of the input.json file
```bash
symetrical_reconversion: bool
ref_name: string or list[string]            # file name of peak list from reference version 
trans_name: string or list[string]           # file name of peak list from transfer version 
dimension: integer                      # dimensionality of spectra
TC: float                           # time of ccr evolution 
NS: list[integer]                       # number of scans in the reference and the transfer version

angle_pos: list[integer]                # angle position of all peaks: -2, -1, 0, 1, 2
CCR_pos:                            # number of measured angle: 1, 2
nucl: list[string]                     # nuclei of all peaks: H, N, C, CA, CB, HA, HB
pos_nucl: list[integer]                 # nuclei position of all peaks: -2, -1, 0, 1, 2 
angle_num: integer                      # number of measured angle: 1, 2
angle_names: list[string]              # angle: phi, psi
noise: list[integer]                    # noise level in the reference and the transfer version
other: string                          # add this as a readable comment for the script to particular experiment
```
Example if experiment is in dictionary (`CCR_dict.py` file) and was recorded withot symetrical reconversion:
```json
{"CCR_1_100NUS":{
    "symetrical_reconversion": false,
    "type_of_CCR": "CCR_1",
    "ref_name": "CCR_1_100NUS_a",
    "trans_name": "CCR_1_100NUS_x",
    "dimension": 4,
    "NS": [4,24],
    "TC": 0.0286,
    "other": "100NUS"
    }
}
```

Example if experiment is NOT in dictionary (`CCR_dict.py` file):
```json
{"CCR_x":{
    "symetrical_reconversion": false,
    "type_of_CCR": "CCR_x",
    "ref_name": "CCR_1_100NUS_a",
    "trans_name": "CCR_1_100NUS_x",
    "dimension": 4,
    "NS": [4,24],
    "TC": 0.028,
    "angle_num": 1,
    "CCR_pos": -1,
    }
}
```
If you are using symetrical reconversion approach you should extend `ref_name`, `trans_name` to list of the transfers and the references versions file names and extend `NS` to four-item list (ref, ref, trans, trans). 
```json
"CO_COCA_trans":{
    "symetrical_reconstrution": true,
    "type_of_CCR": "CO_COCA_trans",
    "ref_name": ["63_CO_COCA_1_1", "66_CO_COCA_2_2"],
    "trans_name": ["64_CO_COCA_1_2","65_CO_COCA_2_1"],
    "dimension": 3,
    "NS": [16,16,32,32],
    "TC": 0.06,
    "angle_num": 1,
    "CCR_pos": -1,
    },
```


### Additionally options:
You can compere different data for the same CCR rate. If the script finds 2 data sets for the same type of CCR rate then it will prepare a graph comparing them. Remember to add comments to specific experiments in the experiment_set file like:
``` 
other: [number - 3 digits][parameter - 3 letters]
``` 
example:  
```
other: 100NUS
```     

### Noise level vs peaks noise

Peak uncertainity is calculated from noise level for each peak. The script can use the same noise level for every peak or use informations from individual noise level for each peak. 
The information about noise level can be place:
- in `input.json` as parametr "noise"
- in peak list with posistions in points - second row "Average noise level = "+number
- in file with peak uncertainity (`name+"_peaks_noise.list"`)


### Output:
There is several types of output:
- table with CCR rates values and standard deviation for every residue for every experiments (`CCRrate.csv` file)
- table for each experiment with CCR rates for every residue (`type_of_CCR+".csv"`)
- text file (`RaportBox.txt`) with whole terminal output and additional informations to evaluate script work



<br>

---  
---

### Citation:  
If you use those scripts please cite:




