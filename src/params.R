params <- 
  list(
    
    starting_ages = 50:74,
    starting_age_distribution = c(
      rep(0.042095795, 5),
      rep(0.045708717, 5),
      rep(0.044164605, 5),
      rep(0.039997002, 5),
      rep(0.028033881, 5)
    ),
    max_screening_age = 75,
    
    proportion_entering_per_year = 1,
    
    low_risk_polyp_model = source('./src/low_risk_polyp_model.R')$value,
    high_risk_polyp_model = source('./src/high_risk_polyp_model.R')$value,
    
    prevalence_preclinical_early_CRC_model = source('./src/preclinical_early_CRC_model.R')$value,
    prevalence_preclinical_regional_CRC_model = source('./src/preclinical_regional_CRC_model.R')$value,
    prevalence_preclinical_distant_CRC_model = source('./src/preclinical_distant_CRC_model.R')$value,
    
    background_mortality_model = source('./src/life_table_to_model.R')$value,
    
    normal_mucosa_to_low_risk_polyp_model = source('./src/normal_mucosa_to_low_risk_polyp_model.R')$value,
    
    low_risk_polyp_to_high_risk_polyp = 0.036,
    high_risk_polyp_to_preclinical_local_CRC = 0.042,
    preclinical_local_CRC_to_preclinical_regional_CRC = 0.17,
    preclinical_regional_CRC_to_preclinical_distant_CRC = 0.10,
    preclinical_local_CRC_to_clinical_local_CRC = 0.17,
    preclinical_regional_CRC_to_clinical_regional_CRC = 0.21,
    preclinical_distant_CRC_to_clinical_distant_CRC = 1,
    
    adherence_FOBT = 0.8,
    adherence_FIT = 0.8,
    adherence_sigmoidoscopy = 0.8,
    adherence_colonoscopy = 0.8,
    adherence_colonoscopy_after_positive_screening = 0.8,
    
    FOBT_sens_low_risk_polyp = 0.10,
    FOBT_sens_high_risk_polyp = 0.24,
    FOBT_sens_CRC = 0.70,
    FOBT_spec = 0.93,
    
    FIT_sens_low_risk_polyp = 0.17,
    FIT_sens_high_risk_polyp = 0.42,
    FIT_sens_CRC = 0.92,
    FIT_spec = 0.90,
    
    colonoscopy_sigmoidoscopy_sens_low_risk_polyp = 0.85,
    colonoscopy_sigmoidoscopy_sens_high_risk_polyp = 0.95,
    colonoscopy_sigmoidoscopy_sens_CRC = 0.95,
    sigmoidoscopy_spec = 0.92,
    prob_neg_sigmoidoscopy_proximal_neoplasm = 0.21,
    
    death_colonoscopy_perforation = 0.012,
    diagnostic_colonoscopy_perforation = 0.0008,
    
    clinical_localized_CRC_5_year_mortality = 0.1,
    clinical_regional_CRC_5_year_mortality = 0.35,
    clinical_distant_CRC_5_year_mortality = 0.92,
    
    QALY_local_regional_cancer = 0.7,
    QALY_distant_cancer = 0.25,
    
    colonoscopy_cost = 40,
    FOBT_cost = 5,
    FIT_cost = 10,
    sigmoidoscopy_cost = 20,
    colonoscopy_perforation_tx_cost = 0,
    local_cancer_tx_cost = 0,
    regional_cancer_tx_cost = 0,
    distant_cancer_tx_cost = 0,
    surveillance_cost = 0,
    
    discount_rate = 0.03
)