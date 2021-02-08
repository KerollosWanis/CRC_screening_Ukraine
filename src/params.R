params <- 
  list(
    
    starting_age = 50,
    
    prevalence_low_risk_polyp_age_50 = 0.2,
    prevalence_low_risk_polyp_age_60 = 0.4,
    prevalence_low_risk_polyp_age_70 = 0.5,
    
    prevalence_high_risk_polyp_age_50 = 0.05,
    prevalence_high_risk_polyp_age_60 = 0.09,
    prevalence_high_risk_polyp_age_70 = 0.16,
    prevalence_high_risk_polyp_age_80 = 0.21,
    
    prevalence_preclinical_early_CRC_age_50 = 0.0024,
    prevalence_preclinical_regional_CRC_age_50 = 0.0012,
    prevalence_preclinical_distant_CRC_age_50 = 0.0004,
    
    background_mortality_model = source('./src/life_table_to_model.R')$value,
    
    normal_mucosa_to_low_risk_polyp_age_50 = 0.00836,
    normal_mucosa_to_low_risk_polyp_age_55 = 0.0099,
    normal_mucosa_to_low_risk_polyp_age_60 = 0.01156,
    normal_mucosa_to_low_risk_polyp_age_65 = 0.0133,
    normal_mucosa_to_low_risk_polyp_age_70 = 0.01521,
    
    low_risk_polyp_to_high_risk_polyp = 0.036,
    high_risk_polyp_to_preclinical_local_CRC = 0.042,
    preclinical_local_CRC_to_preclinical_regional_CRC = 0.17,
    preclinical_regional_CRC_to_preclinical_distant_CRC = 0.10,
    preclinical_local_CRC_to_clinical_local_CRC = 0.17,
    preclinical_regional_CRC_to_clinical_regional_CRC = 0.21,
    preclinical_distant_CRC_to_clinical_distant_CRC = 1,
    
    adherence_FOBT = 0.75,
    adherence_FIT = 0.75,
    adherence_sigmoidoscopy = 0.75,
    adherence_colonoscopy = 0.8,
    adherence_colonoscopy_after_positive_screening = 0.84,
    
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
    
    colonoscopy_cost = 100,
    FOBT_cost = 8,
    FIT_cost = 22,
    sigmoidoscopy_cost = 20,
    colonoscopy_perforation_tx_cost = 500,
    local_cancer_tx_cost = 500,
    regional_cancer_tx_cost = 9000,
    distant_cancer_tx_cost = 20000,
    surveillance_cost = 200,
    
    discount_rate = 0.03
)