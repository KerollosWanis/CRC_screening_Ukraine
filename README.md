# A decision analytic model for CRC screening in low/middle income countries
# Introduction
Here we provide the code to reproduce the analysis described in: 

### Citation

> Welten VM, Wanis KN, Semeniv S, Shabat G, Dabekaussen KF, Davids JS, Beznosenko A, Suprun U, Soeteman DI, Melnitchouk N. Colonoscopy Needs for Implementation of a Colorectal Cancer Screening Program in Ukraine. World Journal of Surgery. 2022 Jul 14:1-1

# Abstract
 


# Organization
- `main.R` — R file which runs the Monte Carlo simulation. Calls on scripts contained in `src`. Expects a subfolder named `output`.
- `src`  — Folder containing scripts including main functions for monte carlo simulation, as well as `params.R` which contains the values of the parameters used in the base case. `life_table_to_model.R` contains life table mortality probabilities for the base case. `model_functions - extended.R` contains model functions which add diagnostic and surveillance colonoscopies to the output.
- `summarize_output.R` — R file which uses the simulated data to produce cost-effectiveness results. 

# Correspondence
If you have any questions, comments, or discover an error, please contact Kerollos Wanis at knwanis@g.harvard.edu or knwanis@gmail.com
