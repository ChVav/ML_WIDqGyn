
The purpose is to estimate the risk for endometrial (WID-qEC) or cervical cancer (WID-qCIN) based on the methylation status of DNA recovered from cervical smear samples, assayed with MethyLight.

The WIDqGYN test is based on 6 different targets to predict CIN and/or EC, 2 additional targets to predict IC or epithelial cell content of a sample, and COL2A1 as the ref gene. 
First publication: https://doi.org/10.1200/JCO.22.00266

This pipeline is taylored to use with Quantstudio 7 and accompanying software.

Please not that starting from v1.2.1, the pipeline assumes following gblock standard curve: STD1 10^5, STD2 10^4, STD3 10^3, STD4 10^2 copy numbers/5uL. <br>
For the old standard curve (shifted towards higher concentrations, i.e. STD1 10^6), please use v1.2.0 of the pipeline for correct back calculation of DNA amount input.

A spin off is the WIDqEC test that will be commercialized and is based on 3 targets of the WIDqGYN to predict CIN and/or EC, and COL2A1 as the ref gene. (https://doi.org/10.1002/ijc.34275)
A separate repo for the Shiny app CalculateWIDqEC used for initial commercialization testing (reproducibility between different labs) is kept, that is tailored to the Quantstudio 5 and accompanying software.

## Full Analysis workflow WID-qGYN in R

### Step 1
GeneratePlateLayout Apps (Shiny only) will create the plate layout for a (full) qWID-GYN test from a sample list drawn in Excel (up to 42 samples), for a qWID-EC test or for a qWID-CIN test (up to 90 samples). <br>
Obviously, your plates should be pipetted accordingly. <br>
Make sure R is installed and in your environment path. <br>
Libraries needed: shiny, readxl, stringr and dplyr. <br>
For Windows: e.g., Double-click GeneratePlatelayoutGYN.bat (or a short-cut thereof), and follow instructions in the browser. <br>
For Mac/Unix: could use the .sh scripts (not yet tested).

### Step 2
The resulting csv file from step 1 can be used to analyze a qPCR run (with correct layout) from a Quantstudio 7 with the Design & Analysis Software 2.5.0 (Applied Biosystms TM).

### Step 3
The CalculatePMRGyn function (Shiny or R) is used to estimate Percentage Methylated Reference (PMR) values for the individual targets.

Users inexperienced with R should run the app version of Step3 in the pipeline: <br>
Double-click CalculatePMRGyn.bat (or a short-cut thereof), and follow instructions in the browser. Here the The WID-qEC (+diagnosis) and WID-qCIN is also calculated, depending on the targets included. <br>

In R, additional functions are implemented so that external standards or a fixed intercept/slow can be used in case the standards on the plate have failed.




