# COVID_ILC: Reproducing the Main Figures

This repository contains analysis scripts for reproducing the main figures in our manuscript, **“Single cell multi-omic analysis identifies key genes differentially expressed in innate lymphoid cells from COVID-19 patients.”**

The scripts use processed study data, including manually gated ILC populations, and the supplied pySCENIC outputs.

**[Manuscript](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1374828)**

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/abhinav2kaushik/COVID_ILC.git
cd COVID_ILC
```

### 2. Download and extract the processed data

Download these two archives from **DOI: XXXXXXXXX**:

- `DataFiles.zip`
- `PyScenic_output.zip`

Place both archives in the main `COVID_ILC` folder and extract them:

```bash
unzip DataFiles.zip
unzip PyScenic_output.zip
```

Ensure that the resulting `DataFiles/` and `PyScenic_output/` directories are located directly inside the main `COVID_ILC` folder.

## Expected Output

Running the scripts will create output folders for the corresponding main figures (e.g., `Figure_1/`). Each folder will contain the generated plots and associated p-values, where applicable, as reported in the study.

### 3. Install the required packages

Before running the scripts, install all packages listed in the supplied session information (`sessionInfo`). Use the recorded package versions to match the original analysis environment as closely as possible.

### 4. Run the analysis scripts

Set the main `COVID_ILC` folder as your working directory and run the scripts corresponding to the figures you wish to reproduce. Keep the extracted directory structure intact so that the scripts can locate their input files.

### Optional: Running pySCENIC

Running pySCENIC (https://pyscenic.readthedocs.io/en/latest/) is optional because the essential input files and precomputed output files are already provided in `PyScenic_output.zip`. These outputs can be used directly for the downstream analyses and figure generation. 
The Python script used to generate the pySCENIC outputs is also included for reference and for users who wish to rerun this step.Please refer pySCENIC.py for script we used to generate PyScenic output. 


## Citation

Kaushik A, Chang I, Han X, He Z, Komlosi ZI, Ji X, Cao S, Akdis CA, Boyd S, Pulendran B, Maecker HT, Davis MM, Chinthrajah RS, DeKruyff RH and Nadeau KC (2024) Single cell multi-omic analysis identifies key genes differentially expressed in innate lymphoid cells from COVID-19 patients. *Front. Immunol.* 15:1374828. doi: [10.3389/fimmu.2024.1374828](https://doi.org/10.3389/fimmu.2024.1374828)

## Contact

**Abhinav Kaushik**  
[akaushik@hsph.harvard.edu](mailto:akaushik@hsph.harvard.edu)
