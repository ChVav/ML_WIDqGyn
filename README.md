
The purpose is to estimate the risk for endometrial (WID-qEC) or cervical cancer (WID-qCIN) based on the methylation status of DNA recovered from cervical smear samples, assayed with MethyLight.

The WIDqGYN test is based on 6 different targets to predict CIN and/or EC, 2 additional targets to predict IC or epithelial cell content of a sample, and COL2A1 as the ref gene. 
First publication: https://doi.org/10.1200/JCO.22.00266

This pipeline is taylored to use with Quantstudio 7 and accompanying software.

Please not that v1.2.1 assumes following gblock standard curve: STD1 10^5, STD2 10^4, STD3 10^3, STD4 10^2 copy numbers/5uL. <br>
For the old standard curve (shifted towards higher concentrations, i.e. STD1 10^6), please use v1.2.0 of the pipeline for correct back calculation of DNA amount input.

A spin off is the WIDqEC test that will be commercialized and is based on 3 targets of the WIDqGYN to predict CIN and/or EC, and COL2A1 as the ref gene. (https://doi.org/10.1002/ijc.34275)
A separate repo for the Shiny app CalculateWIDqEC used for initial commercialization testing (reproducibility between different labs) is kept, that is tailored to the Quantstudio 5 and accompanying software.

## Full Analysis workflow WID-qGYN in R

### Step 1
GeneratePlateLayout function will create the plate layout for a qWID-GYN test from a sample list drawn in Excel (up to 42 samples).

### Step 2
The resulting csv file from step 1 can be used to analyze a qPCR run (with correct layout) from a Quantstudio 7 with the Design & Analysis Software 2.5.0 (Applied Biosystms TM).

### Step 3
Use CalculatePMRGyn function is used to estimate Percentage Methylated Reference (PMR) values for the individual targets.


## Instructions for using the Shiny app
Users inexperienced with R may run the app version of the pipeline.

Make sure R is installed and in your environment path.
Libraries needed: shiny, readxl, stringr and dplyr.

For Step1: Double-click GeneratePlatelayoutGYN.bat (or a short-cut thereof), and follow instructions in the browser. <br>
For Step2: Double-click CalculatePMRGyn.bat (or a short-cut thereof), and follow instructions in the browser. Here the The WID-qEC (+diagnosis) and WID-qCIN is also calculated.

In contrast to the original pipeline, there is no way to use external standards or a fixed intercept/slow in case the standards on the plate have failed.