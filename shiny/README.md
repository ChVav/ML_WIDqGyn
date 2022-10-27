
Make sure R is installed and in your environment path.

Libraries needed: shiny, readxl, stringr and dplyr.

### Step 1
GeneratePlateLayoutGYN will create the plate layout for a qWID-GYN test from a sample list drawn in Excel (up to 42 samples).
Layout version where reps are pipeted underneath each other
Double-click GeneratePlatelayoutGYN.bat (or a short-cut thereof), and follow instructions in the browser.

### Step 2
The resulting csv file from step 1 can be used to analyze a qPCR run (with correct layout) from a Quantstudio 7 with the Design & Analysis Software 2.5.0 (Applied Biosystms TM).

### Step 3
Use CalculatePMRGyn to estimate Percentage Methylated Reference (PMR) values for the individual targets.
Double-click CalculatePMRGyn.bat (or a short-cut thereof), and follow instructions in the browser.