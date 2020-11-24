########################################
# Load dependencies
########################################

source('./src/dependencies.R')

########################################
# Initialize model parameters
########################################

source('./src/params.R')

########################################
# Load model functions
########################################

source('./src/model_functions.R')

########################################
# Set monte carlo simulation values
########################################

N <- 10000000
N_years <- 50
screening_type <- 'none' #choose 'FOBT', 'FOBT + flex sig', 'colonoscopy', or 'none'
all_years <- 'n' #y for list including sim at every year, n for final year only

########################################
# Initialize monte carlo simulation data
########################################

sim_data <- initialize_sim_data(params, N)

########################################
# Run monte carlo simulation
########################################

if (all_years == 'y'){
  sim_data_list <- list()
}

pb <- startpb(min=0, max=N_years)

for(year in 1:N_years) {
  
  sim_data <- update_sim(params, sim_data, screening_type)
  
  if (all_years == 'y'){
    sim_data_list[[year]] <- sim_data
  }
  
  setpb(pb, year)
}

if (all_years == 'y'){
  save(sim_data_list, file=paste0("./output/", screening_type, "_", N, "_", N_years, "_all_years.RData"))
} else {
  save(sim_data, file=paste0("./output/", screening_type, "_", N, "_", N_years, ".RData"))
}