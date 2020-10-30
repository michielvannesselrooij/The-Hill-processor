# The Hill processor
Dimple Aerospace B.V. - MATLAB package for processing and analyzing Hill measurements

## Getting started
1. Add `Processing` folder to your MATLAB path
2. Rename `template_mastersheet.xlsx` to `mastersheet.xlsx` to start using it

## Usage
**Add new measurement data**
  - Inside the `Measurements` folder, create a folder that contains the Hill output files for one set of measurements (reference and one target plate)
  - Include three files per velocity sweep: `.csv` + `_F.csv` + `_p.csv`
  - The alphabetical order of the files determines the order of interpretation
    - *set #1*: ref (if `warmup` set to `yes` in mastersheet file)
    - *set #2*: ref
    - *set #3*: target
    - *set #4*: ref
    - *set #5*: target
    - *set #6*: ref
    - *set #7*: etc...
  - If hot-wire data is available, add it to a folder `hotwire` inside the  folder with csv files
- Add a row to `mastersheet.xlsx` for the measurement. Reference the name of the folder you just created to store the measurement files. The id of this row can be referenced in any analysis to include this measurement in a comparison.

**Create new analysis**
- Inside the `Analyses` folder, duplicate the `_templateAnalysis` and give it a suitable name
- Inside the new folder, open `config.xlsx` and enter the measurement id's from `mastersheet.xlsx`
- *(Optional)* Specify plot styling for each measurement included in the analysis
- Run `analyse.m` in MATLAB. If analysis includes new unprocessed measurements this will take a few minutes. After the first processing, the raw csv measurement files can be removed and archived to save storage space
- Results will be stored as `.png` and `.fig` figures inside the analysis folder

**Testing**
Run the demo analysis `analyse.m` inside `Analyses/_templateAnalysis` and check for errors