
#testing before committing to v1.0.0

source("../../../rpip/pmr.R")

# old thresholds, on plate
test <- pmr(folder="../../../data/dummy/1-input/", 
                output="test1/", 
                "dummytest",
                write.results = TRUE,
                threshold_COL2A1=35,
                threshold_targets=50,
                external_curve=NULL,
                calib_fixed = FALSE,
                fix_intercept = 36.9,
                fix_slope = -3.4,
                path = "defined",
                autotest= FALSE,
                path_lib="C:/Users/charl/Desktop/Charlotte_DevOp/ML_WIDqGYn/lib/")

# read in expected data
data_expected <- read.xlsx("../../../data/dummy/2-output/COL2A1_t35/checkOnplateCalib.xlsx", sheet=1)

# compare results, sanity checks
identical(test[[2]]$Sample,data_expected$Sample) #TRUE
identical(round(test[[2]]$ImC,3),round(data_expected$ImC,3)) #TRUE
identical(round(test[[2]]$DPP6,3),round(data_expected$DPP6,3))
identical(round(test[[2]]$GYPC1,3),round(data_expected$GYPC1,3))

test <- pmr(folder="../../../data/dummy/1-input/", 
            output="test2/", 
            "dummytest",
            write.results = TRUE,
            threshold_COL2A1=35,
            threshold_targets=50,
            external_curve=NULL,
            calib_fixed = TRUE,
            fix_intercept = 36.9,
            fix_slope = -3.4,
            path = "defined",
            autotest= FALSE,
            path_lib="C:/Users/charl/Desktop/Charlotte_DevOp/ML_WIDqGYn/lib/")

# read in expected data
data_expected <- read.xlsx("../../../data/dummy/2-output/COL2A1_t35/checkFixedCoeff.xlsx", sheet=1)

# compare results, sanity checks
identical(test[[2]]$Sample,data_expected$Sample) #TRUE
identical(round(test[[2]]$ImC,3),round(data_expected$ImC,3)) #TRUE
identical(round(test[[2]]$DPP6,3),round(data_expected$DPP6,3))
identical(round(test[[2]]$GYPC1,3),round(data_expected$GYPC1,3))
