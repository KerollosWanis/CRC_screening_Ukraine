# A decision analytic model for CRC screening in low/middle income countries
# Introduction
Here we provide the code to reproduce the analysis described in: 

### Citation

> 

# Abstract
 


# Organization
- `main.R` — R file which runs the Monte Carlo simulation. Calls on scripts contained in `src`. Expects a subfolder named `output`.
- `src`  — Folder containing scripts including main functions for monte carlo simulation, as well as `params.R` which contains the values of the parameters used in the base case. `life_table_to_model.R` contains life table mortality probabilities for the base case. 
- `summarize_output.R` — R file which uses the simulated data to produce cost-effectiveness results. 

# Correspondence
If you have any questions, comments, or discover an error, please contact Kerollos Wanis at knwanis@g.harvard.edu or knwanis@gmail.com
