# GeneratePlateLayoutGYN v0.0.2

Layout version where reps are pipeted underneath each other

## Full Analysis workflow WID-qGYN
The purpose of the WID-qGYN is to estimate the risk for endometrial (WID-qEC) or cervical cancer (WID-qCIN) based on the methylation status of DNA recovered from cervical smear samples, assayed with MethyLight.

### Step 1
GeneratePlateLayoutGYN will create the plate layout for a qWID-GYN test from a sample list drawn in Excel (up to 42 samples).

### Step 2
The resulting csv file from step 1 can be used to analyze a qPCR run (with correct layout) from a Quantstudio 7 with the Design & Analysis Software 2.5.0 (Applied Biosystms TM).

### Step 3
Use CalculatePMR (v.0.0.2) to estimate Percentage Methylated Reference (PMR) values for the individual targets.

## Instructions for GeneratePlateLayoutGYN (Step 1)

Make sure R is installed and in your environment path.

Libraries needed: shiny, readxl, stringr and dplyr.

Double-click GeneratePlatelayoutGYN.bat (or a short-cut thereof), and follow instructions in the browser.