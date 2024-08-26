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
```bash
python3 read_ucsf [-h] [-np NUMBER_OF_POINTS_FOR_NOISE] [-pl PEAK_LEVEL] [-nrm] [-op] [-o OUTPUT_NAME] [-n] [-sn SIGNALTONOISE] [-rec] spectrum_path peak_list_path
```

### Positional arguments:        
```bash
  spectrum_path         path to spectrum file (Sparky format, *.ucsf, https://nmrfam.wisc.edu/nmrfam-sparky-distribution/)
  peak_list_path        path to peak list file (Sparky format)
```

### Optional arguments:      
```bash
  -h, --help            show the help message and exit
  -np, --npoints        the number of randomly chosen points used to calculate the noise level (default: 1000 for 2D spectrum, 10000 for 3D, 100000 for 4D)
  -pl, --peaklevel      the minimal peak height, in scientific notation e.g. 1e+7 (default: XXXXXXXXXXXXXXXXXXXXXXX)
  -nrm, --noRemove      do NOT remove invisible peaks
  -op, --onlypoints     change ppm value to spectral points value
  -o, --output_name     output name (default: XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx)
  -n, --noise           calculate only the noise level
  -sn, --signal3noise   the minimal signal to noise ratio (default: none)
  -rec, --reconstructedspectrum      calculate the spectral noise using the random points at peak proton position - required for spectra reconstructed from NUS data (if this argument is not used, the noise is calculated using random points from across the spectrum)
```
### Input files:
  - multidimensional NMR spectra in ucsf-[Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format
  - peak list in [Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format  

### Output:
The directory with the following files:
- the original peak list with peak positions given in  spectral points (`name+"_origin_points.list"`)
- two peak lists after adjusting peak positions, one contains peak positions in ppm (`name+"_new_ppm.list"`), another one in spectral points (`name+"_new_points.list"`)
- the peak list which contains only peak names and respective noise levels (`name+"_peaks_noise.list"`)
- the text file (`info.txt`) with whole terminal output and additional information to evaluate script performance



<br><br>

## calc_CCR_rate.py

The script calculates the cross-correlated relaxation (CCR) rates, using peak intensities from the two spectra, reference and transfer (quantitative gamma approach). 
   
Our scripts can work with any type of CCR rates but a few of them were already pre-programmed:
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
*Additionally for CCR 5 and 6 the CCR rates are multiplied by -1, to achieve proper results.*


### Input files (obligatory):  
- peak lists with peak positions in ppm and peak heights (`name+"_ppm.list"`)
- peak lists with peak positions in spectral points and peak heights (`name+"_points.list"`)
- JSON file with experiments description ('XXXXXXXXXXXXXXXXXXXXXXXXXXXX')   
- sequence in FASTA format (`seq` file)     

### Input files (optional): 
- peak lists with peak names and noise level (`name+"_peaks_noise.list"`)
- file with reference values of CCR rates (XXXXXXXXXXXXXXXXXXXXXX) - for results validation

<br>

### Launch in command line:    
```bash
python3 calc_CCR_rate [-h] [-s SEQ_FILE_NAME] [-r REFGAMMA] [-e EXPSET] [-pub] [-pres] file_directory
```

### Positional arguments:       
```bash
file_directory         path to the directory with all required files
```

### Optional arguments:       
```bash
-h, --help              show the help message and exit
-s, --seq               name of the file with amino acid sequence (default: 'seq')
-r, --refgamma          name of the file with reference values of CCR rates (file must be in csv format, columns names should be: AA, CCR_name_1, CCR_name_2, CCR_name_2, ...)
-e, --expset            experiments setup file ('input.json') - the structure of the file must be the same as the original one
-pub, --publication     prepare output figures in publication size
-pres, --presentation   prepare output figures in presentation size
```



### Structure of the 'input.json' file

```bash
symmetrical_reconversion: bool          # true: symmetrical reconversion approach
ref_name: string or list[string]        # file name of the peak list of the reference spectrum 
trans_name: string or list[string]      # file name of the peak list of the transfer spectrum 
dimension: integer                      # dimensionality of the spectra
TC: float                               # length of the ccr evolution delay
NS: list[integer]                       # number of scans in the reference and the transfer version of the experiment

angle_pos: list[integer]                # relative angle position: -2, -1, 0, 1, 2
CCR_pos:                                # relative CCR rate position: -2, -1, 0, 1, 2 XXXXXXXXXXXXXXXXXXXXXX nie rozumiem, co ten parametr określa
nucl: list[string]                      # nuclei of all peaks: H, N, C, CA, CB, HA, HB XXXXXXXXXXXXXXXXXXXXX do czego jest ta informacja? W przykładach poniżej tego parametru nie ma...
pos_nucl: list[integer]                 # nuclei position of all peaks: -2, -1, 0, 1, 2 XXXXXXXXXXXXXXXXXXXXXXX j.w.
angle_num: integer                      # number of measured angles: 1, 2
angle_names: list[string]               # angles: phi, psi
noise: list[integer]                    # noise level in the reference and the transfer spectrum
other: string                           # a readable comment for the script to particular experiment
```
Example if an experiment is in the dictionary `CCR_dict.py` and was recorded without symmetrical reconversion: XXXXXXXXXXXXXXXXX nie rozumiem, co to za katalog CCR_dict.py i dlaczego ma rozszerzenie .py? 
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
You can compare different data for the same CCR rate. If the script finds 2 data sets for the same type of CCR rate then it will prepare a graph comparing them. XXXXXXXXXXXXXXXXXX sam to zrobi czy trzeba mu dać plik z 'reference CCR values' (jeden z opcjonalnych) Remember to add comments to specific experiments in the `input.json` like:
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
There are several input files:
- table with CCR rates values and their uncertainties for each residue for each experiment (`CCRrate.csv`) XXXXXXXXXXXXXXX nie rozumiem różnicy pomiędzy tym i kolejnym
- table for each experiment with CCR rates for each residue (`type_of_CCR+".csv"`)
- text file (`RaportBox.txt`) with whole terminal output and additional information to evaluate script performance



<br>

---  
---

### Citation:  
If you use these scripts please cite:




