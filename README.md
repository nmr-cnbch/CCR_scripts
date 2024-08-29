# CCR_scripts   
A set of Python scripts for analyzing NMR spectra  for cross-correlated relaxation (CCR) effect measurements. The set includes:
  - read_ucsf.py (adjusting the positions of peaks in the peak list, to fit the spectral peaks' positions)
  - calc_CCR_rate.py (calculating the CCR rates, using peak intensities)


### software required:
  - [Python 3.6](https://www.python.org/) with libraries:
    - [NumPy](https://www.numpy.org/) 
    - [Matplotlib](https://matplotlib.org/)
    - [SciPy](https://www.scipy.org/)



--- 
---  


## read_ucsf.py    
The script adjusts the positions of peaks in the peak list (provided as an input) using the spectrum file (also provided as input, in [Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format, *.ucsf). The output is the corrected peak list. The peak positions are given in ppm and spectral points.

### Launch in a command line:  
```
python3 read_ucsf.py [-h] [-np NUMBER_OF_POINTS_FOR_NOISE] [-pl PEAK_LEVEL] [-nrm] [-o OUTPUT_DIR] [-sn SIGNALTONOISE] [-rec] spectrum_path peak_list_path
```

### Obligatory arguments:        
```
  spectrum_path         path to spectrum file (Sparky format, *.ucsf, https://nmrfam.wisc.edu/nmrfam-sparky-distribution/)
  peak_list_path        path to peak list file (Sparky format)
```

### Optional arguments:      
```
  -h, --help            show the help message and exit
  -np, --npoints        the number of randomly chosen points used to calculate the noise level (default: 1000 for 2D spectrum, 10000 for 3D, 100000 for 4D)
  -pl, --peaklevel      the minimal peak height, in scientific notation e.g. 1e+7 (default: XXXXXXXXXXXXXXXXXXXXXXX)
  -nrm, --noRemove      do NOT remove invisible peaks (use this option always for transfer versions of experiments!)
  -o, --output_dir      name of output directory (default: spectrum file name from spectrum_path)
  -sn, --signal3noise   the minimal signal-to-noise ratio (default: XXXXXXXXXXXXXXXXXXXXXXXXXXXX)
  -NUS, --NUS           for non-uniformly sampled data, calculate the spectral noise using the random points at the peak's proton position (if this argument is not used, for conventionally-sampled data, the noise is calculated using random points from across the spectrum)
```
### Input files:
  - multidimensional NMR spectra in ucsf-[Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format
  - peak list in [Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format  

### Output:
The directory with the following files:
- the text file (`info.txt`) with the whole terminal output and additional information to evaluate script performance (caution! the content is appended to the file, remove the file if you want to create a new one)
- the original peak list with peak positions given in spectral points (`spectrum_name+"_origin_points.list"`)
- two peak lists after adjusting peak positions, one contains peak positions in ppm (`spectrum_name+"_new_ppm.list"`), another one in spectral points (`spectrum_name+"_new_points.list"`)
- the peak list which contains only peak names and respective noise levels (`spectrum_name+"_peaks_noise.list"`) - appears only if the -NUS option is used

`spectrum_name` is the spectrum file name from spectrum_path


<br><br>

## calc_CCR_rate.py

The script calculates the cross-correlated relaxation (CCR) rates, using peak intensities from the two spectra, reference and transfer (quantitative gamma approach). 
   
Our scripts can work with any CCR rates, but a few of them are pre-defined:
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
$\Gamma = \frac{1}{Tc}\cdot arctanh(\frac{I_{trans}}{I_{ref}}\frac{NS_{ref}}{NS_{trans}})$     
where: $I_{trans/ref}$ - peak intensity in the transfer/reference version, $NS_{trans/ref}$ - number of scans in the transfer/reference version     


### Input files (obligatory):  
- peak lists with peak positions in ppm and peak heights (`name+"_ppm.list"`)
- peak lists with peak positions in spectral points and peak heights (`name+"_points.list"`)
- JSON file with experiments description (default: `input.json`)   
- sequence in FASTA format (default: `seq`)     

### Input files (optional): 
- peak lists with peak names and noise level (`name+"_peaks_noise.list"`)

<br>

### Launch in command line:    
```
python3 calc_CCR_rate.py [-h] [-s SEQ_FILE_NAME] [-e EXPSET] file_directory
```

### Obligatory arguments:       
```
file_directory         path to the directory with all required files
```

### Optional arguments:       
```
-h, --help              show the help message and exit
-s, --seq               name of the file with amino acid sequence (default: `seq`)
-e, --expset            experiments setup file (default: `input.json`)
```



### Structure of the 'input.json' file
```
type_of_CCR                             # name of pre-defined CCR experiment (
CCR_pos *:                              # position of residue (relative to the residue of the directly-detected nucleus) which CCR rate provides information on: -2, -1, 0, 1, 2
dimension **: integer                   # dimensionality of the spectra 
TC ***: float                           # length of the CCR relaxation delay in seconds 
ref_name: string or list[string]        # file name of the peak list of the reference spectrum 
trans_name: string or list[string]      # file name of the peak list of the transfer spectrum 
NS: list[integer]                       # number of scans in the reference and the transfer version of the experiment

```
* if type_of_CCR is given, this value is set automatically
** if type_of_CCR is given, this value is set automatically, however, it can be changed here if eg. the experiment was run in a lower-dimensional version
*** for CCR_1, CCR_2, CCR_5 and CCR_6, this value is set automatically: 0.0286
  
Example of 'input.json' file for experiment pre-defined in `CCR_dict.py`:
```json
{"CCR_1_100NUS":{
    "type_of_CCR": "CCR_1",
    "ref_name": "CCR_1_100NUS_a",
    "trans_name": "CCR_1_100NUS_x",
    "NS": [4,24]
    }
}
```

Example of 'input.json' file for experiment not pre-defined in `CCR_dict.py`:
```json
{"CCR_x":{
    "type_of_CCR": "CCR_x", XXXXXXXXXXXXXXXXXXXXXXXX czy to jest obowiązkowy argument???
    "CCR_pos": -1,
    "dimension": 3,
    "TC": 0.028,
    "ref_name": "CCR_x_100NUS_a",
    "trans_name": "CCR_x_100NUS_x",
    "NS": [4,24]
    }
}
```

### CCR rates uncertainties
CCR rate uncertainties are calculated from the spectral noise levels.
For conventionally sampled spectra:
- the script uses the same noise level for every peak of a given spectrum
- the noise level is saved in the peak list with positions in points - second row "Average noise level = "+number
For spectra reconstructed from non-uniformly sampled (NUS) data:
- the script uses information from the individual noise levels for individual peaks. 
- the noise levels are saved in the file with peak uncertainty (`name+"_peaks_noise.list"`)


### Output:
There are several output files:
- summary file (for all experiments), with all CCR rate values and their uncertainties for each residue (`CCRrate.csv`)
- separate file for each experiment with CCR rates for each residue (`type_of_CCR+".csv"`) - in the 'all_outputs' subdirectory
- text file (`RaportBox.txt`) with whole terminal output and additional information to evaluate script performance - in the 'all_outputs' subdirectory

<br>

---  
---

### Citation:  
If you use these scripts please cite:




