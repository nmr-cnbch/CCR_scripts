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
The script compares UCSF file and peak list (in [Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format), and adjusts the positions of the peaks. The output is the corrected peak lists in ppm value and spectra points.

### Launch in command line:  
```bash
python3 read_ucsf [-h] [-np NUMBER_OF_POINTS_FOR_NOISE] [-pl PEAK_LEVEL] [-nrm] [-op] [-o OUTPUT_NAME] [-n] [-sn SIGNALTONOISE] [-rec] ucsf_path peak_list_path
```

### Positional arguments:        
```bash
  ucsf_path             path to UCSF file
  peak_list_path        path to peak list in [Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format
```

### Optional arguments:      
```bash
  -h, --help            show this help message and exit
  -np, --npoints        to change the number of points used to calculate the noise level: N^(spectra dimentionality + 1); normally is N = 10
  -pl, --peaklevel      add this to set up minimal peak height, in scientific notation e.g. 1e+7
  -nrm, --noRemove      add this if you do NOT want to remove invisible peaks
  -op, --onlypoints     add this if you want to only change ppm value to points value
  -o, --output_name     add this if you want specific output name
  -n, --noise           add this if you want to calculate only the noise level
  -sn, --signal3noise   add this if you want to setup minimal signal to noise ratio
  -rec, --reconstructedspectrum
                        add this if your spectrum was reconstructed. For reconstructed spectra we calculate 'noise' by measuring the height of random points on H-cross-section with peaks. Otherwise, for traditionally recorded spectra the noise is calculated by measuring the height of random points from across the spectrum
```
### Files required:
  - multidimensional NMR spectra in ucsf-[Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format
  - peak list in [Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format  

### Output:
The directory with the following files:
- original peak list with chemical shifts converted to points in the spectrum (`name+"_origin_points.list"`)
- two peak lists after moving peaks, one contains peaks positions in chemical shifts (`name+"_new_ppm.list"`), another one with points in the spectrum (`name+"_new_points.list"`)
- peak list which contains only peak names and uncertainties (`name+"_peaks_noise.list"`)
- text file (`info.txt`) with whole terminal output and additional information to evaluate script functionality



<br><br>

## calc_CCR_rate.py

The script to calculate cross-correlated relaxation (CCR) rates, measured by quantitative approach (two spectra: reference and transfer). 
   
Our scripts can work with any type of CCR rate but a few of them were already pre-programmed:
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


An equation used for calculating CCR rates is:   
$\Gamma = \frac{1}{Tc}\cdot arctanh(\frac{I_{trans}}{I_{ref}})$    
*Additionally for CCR 5 and 6 the CCR rates are multiplied by -1, which is consistent with the experiment.*


### Files required:  
- peak lists with peak positions in chemical shifts (`name+"_ppm.list"`) and peak heights
- peak lists with peak positions in points (`name+"_points.list"`) and peak heights
- JSON file with experiments description   
- sequence in FASTA format (`seq` file)     

### Additional files: 
- peak lists with peak names and uncertainties (`name+"_peaks_noise.list"`)
- file with reference values of CCR rates

<br>

### Launch in command line:    
```bash
python3 calc_CCR_rate [-h] [-s SEQ_FILE_NAME] [-r REFGAMMA] [-e EXPSET] [-pub] [-pres] file_directory
```

### Positional arguments:       
```bash
file_directory         path to directory with all required files
```

### Optional arguments:       
```bash
-h, --help              show this help message and exit
-s, --seq               if name of the file with amino acid sequence is not 'seq', add this with the name of file
-r, --refgamma          if you have a file with reference values of CCR rates, add this with the name of the file (file must be .csv, columns names should be: AA, CCR_name_1, CCR_name_2, CCR_name_2, ...)
-e, --expset            if you want to use experiments setup file with different filename than 'input.json' (structure of the file must be the same as the original file)
-pub, --publication     if you want the picture outputs to be in the publication size
-pres, --presentation   if you want the picture outputs to be in the presentation size
```



### Structure of the 'input.json' file

```bash
symmetrical_reconversion: bool
ref_name: string or list[string]        # file name of the peak list from the reference version 
trans_name: string or list[string]      # file name of the peak list from the transfer version 
dimension: integer                      # dimensionality of the spectra
TC: float                               # time of the ccr evolution 
NS: list[integer]                       # number of scans in the reference and the transfer version

angle_pos: list[integer]                # relative angle position: -2, -1, 0, 1, 2
CCR_pos:                                # relative CCR rate position: -2, -1, 0, 1, 2
nucl: list[string]                      # nuclei of all peaks: H, N, C, CA, CB, HA, HB
pos_nucl: list[integer]                 # nuclei position of all peaks: -2, -1, 0, 1, 2 
angle_num: integer                      # number of measured angles: 1, 2
angle_names: list[string]               # angles: phi, psi
noise: list[integer]                    # noise level in the reference and the transfer version
other: string                           # add this as a readable comment for the script to particular experiment
```
Example if an experiment is in the dictionary `CCR_dict.py` and was recorded without symmetrical reconversion:
```json
{"CCR_1_100NUS":{
    "symmetrical_reconversion": false,
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

If an experiment is NOT in the dictionary `CCR_dict.py`, you have to add information about the position of the CCR rate raletive to directed measure nucleus:
```json
{"CCR_x":{
    "symmetrical_reconversion": false,
    "type_of_CCR": "CCR_x",
    "ref_name": "CCR_1_100NUS_a",
    "trans_name": "CCR_1_100NUS_x",
    "dimension": 4,
    "NS": [4,24],
    "TC": 0.028,
    "CCR_pos": -1,
    }
}
```
If you are using symmetrical reconversion approach you should extend `ref_name`, `trans_name` to a list of the transfer and the reference versions file names and extend `NS` to four-item list (ref, ref, trans, trans). 
```json
"CO_COCA_trans":{
    "symmetrical_reconstrution": true,
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


### Additional options:
You can compare different data for the same CCR rate. If the script finds 2 data sets for the same type of CCR rate then it will prepare a graph comparing them. Remember to add comments to specific experiments in the `input.json` like:
``` 
other: [number - 3 digits][parameter - 3 letters]
``` 
example:  
```
other: 100NUS
```     

### Noise level vs peaks noise

Peak uncertainty is calculated from the noise level for each peak. The script can use the same noise level for every peak or use information from the individual noise level for each peak. 
The information about noise level can be placed:
- in the `input.json` file as the parameter "noise"
- in the peak list with positions in points - second row "Average noise level = "+number
- in the file with peak uncertainty (`name+"_peaks_noise.list"`)


### Output:
There are several types of output:
- table with CCR rates values and standard deviation for every residue for every experiment (`CCRrate.csv`)
- table for each experiment with CCR rates for every residue (`type_of_CCR+".csv"`)
- text file (`RaportBox.txt`) with whole terminal output and additional information to evaluate script functionality



<br>

---  
---

### Citation:  
If you use those scripts please cite:




