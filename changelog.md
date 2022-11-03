
## [1.2.1] - 03.11.2022

## Added
* output additional list checking whether SD for COL2A1 between reps > 1.5 CT, gives warning in QC for WID-qEC
* plot with mean+SD CT COL2A1 for all samples on plate

## [1.1.1] - 27.10.2022 - Released

### Added
* in /rpip, implemented option for using external calibration curve
* this version was to recalibrate the CT thresholds using the Barcelon selfsamples

## [1.0.0] - 26.10.2022

rpip and shiny are now compatible and use the same default thresholds
rpip additionally can be used with precalculated regression coefficients, hence is more flexible

### Fixed
* in /rpip; R pipeline uses correct work functions from /lib

### Deprecated
* in /rpip autotest function does not work anymore, will add this in the future

## [0.2.0] - 26.10.2022
### Added
* in /lib; created separate functions for calibrating COL2A1 gblocks, PMR calculations and WIDqEC interpretation; used in shiny
* in /shiny; added option for CT threshold targets (default 38); COL2A1 default threshold now 30

## [0.1.4] - 25.10.2022
Implementation closest to R pipeline used in https://doi.org/10.1200/JCO.22.00266 and https://doi.org/10.1002/ijc.34275 by Chiara Herzog is in /rpip
* WIDqEC < 0.03 ~ c("low risk EC/CIN", "Negative"); 0.03 =< WIDqEC <= 0.63 ~ c("high risk EC/CIN","Positive"); WIDqEC > 0.63 ~ c("very high risk EC/CIN","Positive")

### Changed
* COL2A1 Ct threshold=35 (36 or higher is kicked out)

### Fixed
In /shiny only
* Using multipoint regression coefficients from COL2A1 standard
* When COL2A1 fails in two reps, PMR is NA, not 0

### Added
* Outputs Excel file batch results: PMR + WID-qEC and WID-qCIN outcome (correctly rounded and formatted), Cq mean, Cq SD, DNA conc estimate, low input (should have negative controls)
* Reprocessing needed: samples where in both reps COL2A1 failed
* Reprocessing recommended: samples where targets did not amplify consistently between reps

