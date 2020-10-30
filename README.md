# The Hill processor
Dimple Aerospace B.V.

## Getting started
1. Add *Processing* folder to your MATLAB path
2. Rename *template_mastersheet.xlsx* to *mastersheet.xlsx* to start using it

## Usage
**Store measurement data**
  - Inside the *Measurements* folder, create a folder that contains the Hill output files for one set of measurements (reference and one target plate)
  - Include three files per sweep: *.csv* + *_F.csv* + *_p.csv*
  - The alphabetical order of the files determines the order of interpretation
    - 1: ref (if *warmup* set to *yes* in mastersheet file)
    - 2: ref
    - 3: target
    - 4: ref
    - 5: target
    - 6: ref
    - 7: etc...
  - If hot-wire data is available, add it to a folder *hotwire* inside the  folder with csv files
- Add a row to *mastersheet.xlsx* for the measurement. Reference the name of the folder you just created to store the measurement files

**Create new analysis**
- Inside the *Analyses* folder, duplicate the *_templateAnalysis* and give it a suitable name
- In 'config.xlsx', enter the measurement id's from *mastersheet.xlsx*