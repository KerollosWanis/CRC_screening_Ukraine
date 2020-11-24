initialize_sim_data <- function(params, N){
  
  #add sample size to params list
  params$N <- N
  
  with(params,{
    
    out <- data.frame(
      #simulation id 
      sim_id = 1:N,
      #starting age
      age = rep(50, N) %>% as.integer(),
      #clinical CRC status variables set to 0 for starting cohort
      clinical_localized_CRC = rep(0, N) %>% as.integer(),
      clinical_regional_CRC = rep(0, N) %>% as.integer(),
      clinical_distant_CRC = rep(0, N) %>% as.integer(),
      clinical_localized_CRC_trt = rep(0, N) %>% as.integer(),
      clinical_regional_CRC_trt = rep(0, N) %>% as.integer(),
      clinical_distant_CRC_trt = rep(0, N) %>% as.integer(),
      surveillance = rep(0, N) %>% as.integer(),
      years_since_dx = rep(0, N) %>% as.integer(),
      dead = rep(0, N) %>% as.integer(),
      CRC_death = rep(0, N) %>% as.integer(),
      #indicators to track screening procedures
      colonoscopy = rep(0, N) %>% as.integer(),
      perforation = rep(0, N) %>% as.integer(),
      FOBT = rep(0, N) %>% as.integer(),
      FOBT_positive = rep(0, N) %>% as.integer(),
      sigmoidoscopy = rep(0, N) %>% as.integer(),
      sigmoidoscopy_positive = rep(0, N) %>% as.integer(),
      #QALYs gained since age 50
      QALYs_gained = rep(0, N) %>% as.integer(),
      #cost incurred during clinical course
      cost_attained = rep(0, N) %>% as.integer(),
      #starting prevalence of low risk polyp
      low_risk_polyp = rbinom(n=N, size=1, prob=prevalence_low_risk_polyp_age_50)
    ) %>% mutate(
      #starting prevalence of other preclinical tumors (probabilities are converted to conditional probabilities)
      high_risk_polyp = case_when(
        low_risk_polyp == 1 ~ as.integer(0),
        TRUE ~ rbinom(n=n(), size=1, 
                      prob=prevalence_high_risk_polyp_age_50/
                        (1-prevalence_low_risk_polyp_age_50))
      ),
      preclinical_early_CRC = case_when(
        low_risk_polyp == 1 | high_risk_polyp == 1 ~ as.integer(0),
        TRUE ~ rbinom(n=n(), size=1, 
                      prob=prevalence_preclinical_early_CRC_age_50/
                        (1-prevalence_low_risk_polyp_age_50-
                           prevalence_high_risk_polyp_age_50))
      ),
      preclinical_regional_CRC = case_when(
        low_risk_polyp == 1 | high_risk_polyp == 1 | preclinical_early_CRC == 1 ~ as.integer(0),
        TRUE ~ rbinom(n=n(), size=1, 
                      prob=prevalence_preclinical_regional_CRC_age_50/
                        (1-prevalence_low_risk_polyp_age_50-
                           prevalence_high_risk_polyp_age_50-
                           prevalence_preclinical_early_CRC_age_50))
      ),
      preclinical_distant_CRC = case_when(
        low_risk_polyp == 1 | high_risk_polyp == 1 | preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 ~ as.integer(0),
        TRUE ~ rbinom(n=n(), size=1, 
                      prob=prevalence_preclinical_distant_CRC_age_50/
                        (1-prevalence_low_risk_polyp_age_50-
                           prevalence_high_risk_polyp_age_50-
                           prevalence_preclinical_early_CRC_age_50-
                           prevalence_preclinical_regional_CRC_age_50))
      ),
      normal_mucosa = case_when(
        low_risk_polyp == 1 | high_risk_polyp == 1 | preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 | preclinical_distant_CRC == 1 ~ as.integer(0),
        TRUE ~ as.integer(1)
      )
    )
    
    return(out)
    
  })
  
}


update_sim <- function(params, sim_data, screening_type) {
  
  #add screening type to params list
  params$screening_type <- screening_type
  
  with(params, {
    
    #update background deaths
    sim_data <- apply_background_death_risk(params, sim_data)
    
    #update QALYs gained
    sim_data <- apply_QALYs_gained(params, sim_data)
    
    #update treatment
    sim_data <- apply_CRC_treatment(params, sim_data)
    
    #update CRC deaths
    sim_data <- apply_clinical_CRC_death_risk(params, sim_data)
    
    #update screening
    sim_data <- apply_screening(params, sim_data)
    
    #update polyp/CRC status due to progression
    sim_data <- apply_CRC_progression(params, sim_data)
    
    #update accumilated cost
    sim_data <- apply_cost(params, sim_data)
    
    #update age and years since dx
    sim_data <- sim_data %>% mutate(
      age = case_when(dead == 1 ~ age,
                      TRUE ~ as.integer(age + 1)),
      years_since_dx = case_when(dead == 0 & 
                                   (clinical_localized_CRC == 1 | clinical_regional_CRC == 1 | clinical_distant_CRC == 1) & surveillance == 1 ~ 
                                   as.integer(years_since_dx + 1),
                                 dead == 1 ~ years_since_dx,
                                 TRUE ~ as.integer(0))
    )
    
    return(sim_data)
    
  })
  
}


apply_background_death_risk <- function(params, sim_data){
  
  with(params, {
    
    #update background deaths
    out <- sim_data %>%
      mutate(
        risk = predict(background_mortality_model, newdata = data.frame(age=age)),
        dead = case_when(
          dead == 1 ~ as.integer(1),
          dead == 0 ~ rbinom(n(), size=1, prob=risk)
        )
        
      ) %>% select(-risk)
    
    return(out)
    
  })
  
}


apply_QALYs_gained <- function(params, sim_data){
  
  with(params, {
    
    #add QALYs based on clinical CRC status
    out <- sim_data %>%
      mutate(
        QALYs_gained = case_when(
          dead == 1 ~ as.numeric(QALYs_gained),
          dead == 0 & (clinical_localized_CRC == 1 | clinical_regional_CRC == 1) ~ QALYs_gained + (QALY_local_regional_cancer*((1-discount_rate)^(age-50))),
          dead == 0 & clinical_distant_CRC == 1 ~ QALYs_gained + (QALY_distant_cancer*((1-discount_rate)^(age-50))),
          TRUE ~ QALYs_gained + (1*((1-discount_rate)^(age-50)))
        )

      )
    
    return(out)
    
  })
  
}


apply_CRC_treatment <- function(params, sim_data){
  
  with(params, {
    
    out <- sim_data %>% mutate(
      
      #apply treatment and indicate transition to surveillance
      clinical_localized_CRC_trt = case_when(
        clinical_localized_CRC == 1 & surveillance == 0 & dead == 0 ~ as.integer(1),
        TRUE ~ as.integer(0)
      ),
      clinical_regional_CRC_trt = case_when(
        clinical_regional_CRC == 1 & surveillance == 0 & dead == 0 ~ as.integer(1),
        TRUE ~ as.integer(0)
      ),
      clinical_distant_CRC_trt = case_when(
        clinical_distant_CRC == 1 & dead == 0 ~ as.integer(1),
        TRUE ~ as.integer(0)
      ),
      surveillance = case_when(
        (surveillance == 1 | clinical_regional_CRC == 1 | clinical_localized_CRC == 1) & dead == 0 & 
          years_since_dx < 5 & 
          (years_since_dx != 0 | clinical_localized_CRC_trt == 1 | clinical_regional_CRC_trt == 1) ~ as.integer(1),
        TRUE ~ as.integer(0)
      )
    )
    
    return(out)
    
  })
  
}


apply_clinical_CRC_death_risk <- function(params, sim_data){
  
  with(params, {
    out <- sim_data %>% mutate(
      #apply risk of death for individuals with clinical CRC
      temp_dead = dead,
      dead = case_when(
        dead == 1 ~ as.integer(1),
        clinical_localized_CRC == 1 ~ rbinom(n(), size=1, prob=1-exp(log(1-clinical_localized_CRC_5_year_mortality)/5)),
        clinical_regional_CRC == 1 ~ rbinom(n(), size=1, prob=1-exp(log(1-clinical_regional_CRC_5_year_mortality)/5)),
        clinical_distant_CRC == 1 ~ rbinom(n(), size=1, prob=1-exp(log(1-clinical_distant_CRC_5_year_mortality)/5)),
        TRUE ~ as.integer(0)
      ),
      CRC_death = case_when(
        CRC_death == 1 ~ as.integer(1),
        CRC_death == 0 & temp_dead == 0 & dead == 1 ~ as.integer(1),
        TRUE ~ as.integer(0)
      )
    ) %>% select(-temp_dead)
    
    return(out)
    
  })
  
}


apply_screening <- function(params, sim_data){
  
  with(params, {

    out <- sim_data %>% mutate(
      
      #evaluate whether patient undergoes FOBT
      FOBT = case_when(screening_type %in% c('FOBT', 'FOBT + flex sig') & dead == 0 & age <= 75 &
                         (normal_mucosa == 1 | low_risk_polyp == 1 | high_risk_polyp == 1 | 
                            preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 | preclinical_distant_CRC == 1) ~ 
                         rbinom(n(), size=1, prob=adherence_FOBT),
                       TRUE ~ as.integer(0)
                       ),
      #evaluate whether patient undergoes flex sig
      sigmoidoscopy = case_when(screening_type == 'FOBT + flex sig' & age%%5 == 0 & dead == 0 & age <= 75 &
                                  (normal_mucosa == 1 | low_risk_polyp == 1 | high_risk_polyp == 1 | 
                                     preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 | preclinical_distant_CRC == 1) ~
                                  rbinom(n(), size=1, prob=adherence_sigmoidoscopy_with_FOBT),
                                TRUE ~ as.integer(0))
    )
    
    out <- out %>% mutate(
      
      #evaluate whether patient has positive test result
      FOBT_positive = case_when(FOBT == 1 & normal_mucosa == 1 ~ rbinom(n(), size=1, prob=1-FOBT_spec),
                                FOBT == 1 & low_risk_polyp == 1 ~ rbinom(n(), size=1, prob=FOBT_sens_low_risk_polyp),
                                FOBT == 1 & high_risk_polyp == 1 ~ rbinom(n(), size=1, prob=FOBT_sens_high_risk_polyp),
                                FOBT == 1 & (preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 | preclinical_distant_CRC == 1) ~ 
                                  rbinom(n(), size=1, prob=FOBT_sens_CRC),
                                TRUE ~ as.integer(0)),
      sigmoidoscopy_positive = case_when(sigmoidoscopy == 1 & normal_mucosa == 1 ~ rbinom(n(), size=1, prob=1-colonoscopy_sigmoidoscopy_spec),
                                         sigmoidoscopy == 1 & low_risk_polyp == 1 ~ rbinom(n(), size=1, prob=colonoscopy_sigmoidoscopy_sens_low_risk_polyp*(1-prob_neg_sigmoidoscopy_proximal_neoplasm)),
                                         sigmoidoscopy == 1 & high_risk_polyp == 1 ~ rbinom(n(), size=1, prob=colonoscopy_sigmoidoscopy_sens_high_risk_polyp*(1-prob_neg_sigmoidoscopy_proximal_neoplasm)),
                                         sigmoidoscopy == 1 & (preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 | preclinical_distant_CRC == 1) ~ 
                                           rbinom(n(), size=1, prob=colonoscopy_sigmoidoscopy_sens_CRC*(1-prob_neg_sigmoidoscopy_proximal_neoplasm)),
                                         TRUE ~ as.integer(0)),
      
    )
    
    #Evaluate whether patient undergoes a colonoscopy based on the screening strategy,
    #the patient's crc status, and the test characteristics
    out <- out %>% mutate(
                              #colonoscopy screening strategy
      colonoscopy = case_when(screening_type == 'colonoscopy' & dead == 0 & age%%10 == 0 & age <= 75 &
                                clinical_localized_CRC == 0 & clinical_regional_CRC == 0 & clinical_distant_CRC == 0 ~ 
                                rbinom(n(), size=1, prob=adherence_colonoscopy),
                              
                              #Other strategies
                              FOBT_positive == 1 | sigmoidoscopy_positive == 1 ~ 
                                rbinom(n(), size=1, prob=adherence_colonoscopy_after_positive_screening), 
                              
                              #others do not get colonoscopy
                              TRUE ~ as.integer(0)
                              )
    )
    
    out <- apply_colonoscopy_risk(params, out)
    
    out <- apply_colonoscopy_benefit(params, out)
    
    return(out)
    
  })
  
}


apply_colonoscopy_risk <- function(params, sim_data){
  
  with(params, {
    
    out <- sim_data %>% mutate(
      perforation = case_when(
        dead == 1 | colonoscopy == 0 ~ as.integer(0),
        colonoscopy == 1 & dead == 0 ~ rbinom(n(), size=1, prob=diagnostic_colonoscopy_perforation),
      )
    )
    
    out <- sim_data %>% mutate(
      dead = case_when(
        perforation == 1 ~ rbinom(n(), size=1, prob=death_colonoscopy_perforation),
        TRUE ~ dead
      )
    )
    
    return(out)
    
  })
  
}


apply_colonoscopy_benefit <- function(params, sim_data){
  
  with(params, {
    
    out <- sim_data %>% mutate(
      
      #remove polyps when undergoing colonoscopy
      normal_mucosa = case_when(
        colonoscopy == 0 | dead == 1 ~ normal_mucosa,
        colonoscopy == 1 & dead == 0 & normal_mucosa == 1 ~ as.integer(1),
        colonoscopy == 1 & dead == 0 & low_risk_polyp == 1 ~ rbinom(n(), size=1, prob=colonoscopy_sigmoidoscopy_sens_low_risk_polyp),
        colonoscopy == 1 & dead == 0 & high_risk_polyp == 1 ~ rbinom(n(), size=1, prob=colonoscopy_sigmoidoscopy_sens_high_risk_polyp),
        colonoscopy == 1 & dead == 0 & (preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 | preclinical_distant_CRC == 1) ~ 
          as.integer(0)
      ),
      
      #clinically identify CRC when screened with preclinical CRC
      clinical_localized_CRC = case_when(
        colonoscopy == 0 | dead == 1 ~ clinical_localized_CRC,
        colonoscopy == 1 & dead == 0 & normal_mucosa == 1 ~ as.integer(0),
        colonoscopy == 1 & dead == 0 & (low_risk_polyp == 1 | high_risk_polyp == 1) ~ as.integer(0),
        colonoscopy == 1 & dead == 0 & preclinical_early_CRC == 1 ~ rbinom(n(), size=1, prob=colonoscopy_sigmoidoscopy_sens_CRC),
        colonoscopy == 1 & dead == 0 & preclinical_regional_CRC == 1 ~ as.integer(0),
        colonoscopy == 1 & dead == 0 & preclinical_distant_CRC == 1 ~ as.integer(0)
      ),
      clinical_regional_CRC = case_when(
        colonoscopy == 0 | dead == 1 ~ clinical_regional_CRC,
        colonoscopy == 1 & dead == 0 & normal_mucosa == 1 ~ as.integer(0),
        colonoscopy == 1 & dead == 0 & (low_risk_polyp == 1 | high_risk_polyp == 1) ~ as.integer(0),
        colonoscopy == 1 & dead == 0 & preclinical_early_CRC == 1 ~ as.integer(0),
        colonoscopy == 1 & dead == 0 & preclinical_regional_CRC == 1 ~ rbinom(n(), size=1, prob=colonoscopy_sigmoidoscopy_sens_CRC),
        colonoscopy == 1 & dead == 0 & preclinical_distant_CRC == 1 ~ as.integer(0)
      ),
      clinical_distant_CRC = case_when(
        colonoscopy == 0 | dead == 1 ~ clinical_distant_CRC,
        colonoscopy == 1 & dead == 0 & normal_mucosa == 1 ~ as.integer(0),
        colonoscopy == 1 & dead == 0 & (low_risk_polyp == 1 | high_risk_polyp == 1) ~ as.integer(0),
        colonoscopy == 1 & dead == 0 & preclinical_early_CRC == 1 ~ as.integer(0),
        colonoscopy == 1 & dead == 0 & preclinical_regional_CRC == 1 ~ as.integer(0),
        colonoscopy == 1 & dead == 0 & preclinical_distant_CRC == 1 ~ rbinom(n(), size=1, prob=colonoscopy_sigmoidoscopy_sens_CRC)
      )
      
    )
    
    
    #update polyp and preclinical crc variables depending on the results of colonoscopy
    out <- out %>% mutate(
      
      low_risk_polyp = case_when(
        normal_mucosa == 1 ~ as.integer(0),
        TRUE ~ low_risk_polyp
      ),
      high_risk_polyp = case_when(
        normal_mucosa == 1 ~ as.integer(0),
        TRUE ~ high_risk_polyp
      ),
      preclinical_early_CRC = case_when(
        clinical_localized_CRC == 1 ~ as.integer(0),
        TRUE ~ preclinical_early_CRC
      ),
      preclinical_regional_CRC = case_when(
        clinical_regional_CRC == 1 ~ as.integer(0),
        TRUE ~ preclinical_regional_CRC
      ),
      preclinical_distant_CRC = case_when(
        clinical_distant_CRC == 1 ~ as.integer(0),
        TRUE ~ preclinical_distant_CRC
      )
      
    )
      
    return(out)
    
  })
  
}


apply_CRC_progression <- function(params, sim_data){
  
  with(params, {
    
    out <- sim_data %>% mutate(
      
      #progress to clinical CRC from preclinical
      clinical_distant_CRC = case_when(
        dead == 0 & preclinical_distant_CRC == 1 ~ rbinom(n(), size=1, prob=preclinical_distant_CRC_to_clinical_distant_CRC),
        TRUE ~ clinical_distant_CRC
      ),
      clinical_regional_CRC = case_when(
        dead == 0 & preclinical_regional_CRC == 1 ~ rbinom(n(), size=1, prob=preclinical_regional_CRC_to_clinical_regional_CRC),
        TRUE ~ clinical_regional_CRC
      ),
      clinical_localized_CRC = case_when(
        dead == 0 & preclinical_early_CRC == 1 ~ rbinom(n(), size=1, prob=preclinical_regional_CRC_to_clinical_regional_CRC),
        TRUE ~ clinical_localized_CRC
      )
    )
    
    out <- out %>% mutate(
      
      #progress into later preclinical stages and update preclinical CRC variables
      #based on results of progression to clinical CRC
      preclinical_distant_CRC = case_when(
        dead == 0 & preclinical_distant_CRC == 1 & clinical_distant_CRC == 1 ~ as.integer(0),
        dead == 0 & preclinical_regional_CRC == 1 & clinical_regional_CRC == 0 ~ 
          rbinom(n(), size=1, prob=preclinical_regional_CRC_to_preclinical_distant_CRC/(1-preclinical_regional_CRC_to_clinical_regional_CRC)),
        TRUE ~ preclinical_distant_CRC
      ),
      preclinical_regional_CRC = case_when(
        dead == 0 & preclinical_regional_CRC == 1 & clinical_regional_CRC == 1 ~ as.integer(0),
        dead == 0 & preclinical_regional_CRC == 1 & preclinical_distant_CRC == 1 ~ as.integer(0),
        dead == 0 & preclinical_early_CRC == 1 & clinical_localized_CRC == 0 ~ 
          rbinom(n(), size=1, prob=preclinical_local_CRC_to_preclinical_regional_CRC/(1-preclinical_local_CRC_to_clinical_local_CRC)),
        TRUE ~ preclinical_regional_CRC
      ),
      preclinical_early_CRC = case_when(
        dead == 0 & preclinical_early_CRC == 1 & clinical_localized_CRC == 1 ~ as.integer(0),
        dead == 0 & preclinical_early_CRC == 1 & preclinical_regional_CRC == 1 ~ as.integer(0),
        dead == 0 & high_risk_polyp == 1 ~ rbinom(n(), size=1, prob=high_risk_polyp_to_preclinical_local_CRC),
        TRUE ~ preclinical_early_CRC
      )
    )
    
    out <- out %>% mutate(
      
      #progress polyps and update polyp variables based on results of preclinical lesion progression
      high_risk_polyp = case_when(
        dead == 0 & high_risk_polyp == 1 & preclinical_early_CRC == 1 ~ as.integer(0),
        dead == 0 & low_risk_polyp == 1 ~ rbinom(n(), size=1, low_risk_polyp_to_high_risk_polyp),
        TRUE ~ high_risk_polyp
      ),
      low_risk_polyp = case_when(
        dead == 0 & low_risk_polyp == 1 & high_risk_polyp == 1 ~ as.integer(0),
        dead == 0 & normal_mucosa == 1 & age >= 50 & age < 55 ~ rbinom(n(), size=1, prob=normal_mucosa_to_low_risk_polyp_age_50),
        dead == 0 & normal_mucosa == 1 & age >= 55 & age < 60 ~ rbinom(n(), size=1, prob=normal_mucosa_to_low_risk_polyp_age_55),
        dead == 0 & normal_mucosa == 1 & age >= 60 & age < 65 ~ rbinom(n(), size=1, prob=normal_mucosa_to_low_risk_polyp_age_60),
        dead == 0 & normal_mucosa == 1 & age >= 65 & age < 70 ~ rbinom(n(), size=1, prob=normal_mucosa_to_low_risk_polyp_age_65),
        dead == 0 & normal_mucosa == 1 & age >= 70 ~ rbinom(n(), size=1, prob=normal_mucosa_to_low_risk_polyp_age_70),
        TRUE ~ low_risk_polyp
      ),
      normal_mucosa = case_when(
        dead == 0 & normal_mucosa == 1 & low_risk_polyp == 1 ~ as.integer(0),
        TRUE ~ normal_mucosa
      )
    )
    
    out <- out %>% mutate(
      
      #return to normal mucosa after 5 years of surveillance in regional or local CRC
      normal_mucosa = case_when(
        dead == 0 & years_since_dx == 5 & (clinical_localized_CRC == 1 | clinical_regional_CRC == 1) ~ as.integer(1),
        TRUE ~ normal_mucosa
      ),
      clinical_regional_CRC = case_when(
        dead == 0 & years_since_dx == 5 & clinical_regional_CRC == 1 ~ as.integer(0),
        TRUE ~ clinical_regional_CRC
      ),
      clinical_localized_CRC = case_when(
        dead == 0 & years_since_dx == 5 & clinical_localized_CRC == 1 ~ as.integer(0),
        TRUE ~ clinical_localized_CRC
      )
    )
      
    return(out)
    
  })
  
}


apply_cost <- function(params, sim_data){
  
  with(params, {
    
    out <- sim_data %>% mutate(
      cost_attained = cost_attained + ((colonoscopy_cost*colonoscopy + FOBT_cost*FOBT +
        sigmoidoscopy_cost*sigmoidoscopy + colonoscopy_perforation_tx_cost*perforation + 
        local_cancer_tx_cost*clinical_localized_CRC_trt + regional_cancer_tx_cost*clinical_regional_CRC_trt +
        distant_cancer_tx_cost*clinical_distant_CRC_trt + surveillance_cost*surveillance)*((1-discount_rate)^(age-50)))
    )
    
    return(out)
    
  })
  
}