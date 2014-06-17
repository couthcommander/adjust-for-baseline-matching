adjust-for-baseline-matching
============================

Role of Matching When Adjusting for Baseline Differences

## R packages to install

    !r
    install.packages(c('lmtest', 'sandwich', 'data.table'))
    # nbpMatching is available on CRAN, but the most latest version is on R-Forge
    install.packages("nbpMatching", repos="http://R-Forge.R-project.org")

## Running the simulation

Edit simulation.R to set an appropriate number of loops.
Line 21 controls outer loops, line 26 controls inner loops.
The defaults (1 and 10) are purposefully small, but will allow the script to run to completion in about 30 minutes.
The values used in our simulation were 100 (outer) and 1000 (inner), which will take closer to two weeks to run to completion.
