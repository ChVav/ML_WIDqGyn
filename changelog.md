
v0.1.4
Implementation closest to R pipeline used in https://doi.org/10.1200/JCO.22.00266 and https://doi.org/10.1002/ijc.34275 is in /rpip
Additions:
* COL2A1 Ct threshold=35 (36 or higher is kicked out)
* WIDqEC < 0.03 ~ c("low risk EC/CIN", "Negative"); 0.03 =< WIDqEC <= 0.63 ~ c("high risk EC/CIN","Positive"); WIDqEC > 0.63 ~ c("very high risk EC/CIN","Positive")

Few fixes/add-ons done in /shiny only
* Using multipoint regression coefficients from COL2A1 standard
* When COL2A1 fails in two reps, PMR is NA, not 0
* Outputs Excel file batch results: PMR + WID-qEC and WID-qCIN outcome (correctly rounded and formatted), Cq mean, Cq SD, DNA conc estimate, low input (should have negative controls)
* Reprocessing needed: samples where in both reps COL2A1 failed
* Reprocessing recommended: samples where targets did not amplify consistently between reps

