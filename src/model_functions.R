initialize_sim_data <- function(params, N){
  
  #add sample size to params list
  params$N <- N
  
  with(params,{
    
    out <- data.frame(
      stringsAsFactors=FALSE,
      #simulation id 
      sim_id = 1:N,
      #starting age
      starting_age = if(length(starting_ages)==1){rep(starting_ages, N)} else {
        sample(starting_ages, size=N, prob=starting_age_distribution, replace=T) %>% as.integer()
      },
      #starting year,
      start_program_year = sample(1:(1/proportion_entering_per_year), size=N, prob=rep(proportion_entering_per_year, 1/proportion_entering_per_year), replace=T) %>% as.integer(),
      #individual screening program type
      individual_screening_type = rep("none", N) %>% as.character(),
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
      #variables to track screening procedures
      hx_of_polyps = rep(0, N) %>% as.integer(),
      hx_of_CRC = rep(0, N) %>% as.integer(),
      colonoscopy = rep(0, N) %>% as.integer(),
      perforation = rep(0, N) %>% as.integer(),
      negative_colonoscopy = rep(0, N) %>% as.integer(),
      low_risk_polyp_identified = rep(0, N) %>% as.integer(),
      high_risk_polyp_identified = rep(0, N) %>% as.integer(),
      FOBT = rep(0, N) %>% as.integer(),
      FIT = rep(0, N) %>% as.integer(),
      FOBT_positive = rep(0, N) %>% as.integer(),
      FIT_positive = rep(0, N) %>% as.integer(),
      sigmoidoscopy = rep(0, N) %>% as.integer(),
      sigmoidoscopy_positive = rep(0, N) %>% as.integer(),
      #QALYs gained since starting age
      QALYs_gained = rep(0, N) %>% as.integer(),
      #cost incurred during clinical course
      cost_attained = rep(0, N) %>% as.integer(),
      #tests performed during clinical course
      num_colonoscopies = rep(0, N) %>% as.integer(),
      num_sigmoidoscopies = rep(0, N) %>% as.integer(),
      num_FIT = rep(0, N) %>% as.integer(),
      num_FOBT = rep(0, N) %>% as.integer()
    ) %>% mutate(
      #variables to track screening procedures
      next_possible_sigmoidoscopy_age = starting_age-1 %>% as.integer(),
      next_possible_colonoscopy_age = starting_age-1 %>% as.integer(),
      age = starting_age,
      #temporary individual prevalence variables
      prevalence_low_risk_polyp = predict(low_risk_polyp_model, newdata = data.frame(age=age)),
      prevalence_high_risk_polyp = predict(high_risk_polyp_model, newdata = data.frame(age=age)),
      prevalence_preclinical_early_CRC = predict(prevalence_preclinical_early_CRC_model, newdata = data.frame(age=age)),
      prevalence_preclinical_regional_CRC = predict(prevalence_preclinical_regional_CRC_model, newdata = data.frame(age=age)),
      prevalence_preclinical_distant_CRC = predict(prevalence_preclinical_distant_CRC_model, newdata = data.frame(age=age)),
      #starting prevalence of low risk polyp
      low_risk_polyp = rbinom(n=n(), size=1, prob=prevalence_low_risk_polyp),
      #starting prevalence of other preclinical tumors (probabilities are converted to conditional probabilities)
      high_risk_polyp = case_when(
        low_risk_polyp == 1 ~ as.integer(0),
        TRUE ~ rbinom(n=n(), size=1, 
                      prob=prevalence_high_risk_polyp/
                        (1-prevalence_low_risk_polyp))
      ),
      preclinical_early_CRC = case_when(
        low_risk_polyp == 1 | high_risk_polyp == 1 ~ as.integer(0),
        TRUE ~ rbinom(n=n(), size=1, 
                      prob=prevalence_preclinical_early_CRC/
                        (1-prevalence_low_risk_polyp-
                           prevalence_high_risk_polyp))
      ),
      preclinical_regional_CRC = case_when(
        low_risk_polyp == 1 | high_risk_polyp == 1 | preclinical_early_CRC == 1 ~ as.integer(0),
        TRUE ~ rbinom(n=n(), size=1, 
                      prob=prevalence_preclinical_regional_CRC/
                        (1-prevalence_low_risk_polyp-
                           prevalence_high_risk_polyp-
                           prevalence_preclinical_early_CRC))
      ),
      preclinical_distant_CRC = case_when(
        low_risk_polyp == 1 | high_risk_polyp == 1 | preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 ~ as.integer(0),
        TRUE ~ rbinom(n=n(), size=1, 
                      prob=prevalence_preclinical_distant_CRC/
                        (1-prevalence_low_risk_polyp-
                           prevalence_high_risk_polyp-
                           prevalence_preclinical_early_CRC-
                           prevalence_preclinical_regional_CRC))
      ),
      normal_mucosa = case_when(
        low_risk_polyp == 1 | high_risk_polyp == 1 | preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 | preclinical_distant_CRC == 1 ~ as.integer(0),
        TRUE ~ as.integer(1)
      )
    ) %>% dplyr::select(-c(prevalence_low_risk_polyp,
                           prevalence_high_risk_polyp,
                           prevalence_preclinical_early_CRC,
                           prevalence_preclinical_regional_CRC,
                           prevalence_preclinical_distant_CRC))
    
    return(out)
    
  })
  
}


update_sim <- function(params, sim_data, screening_type) {
  
  #add screening type to params list
  params$screening_type <- screening_type
  
  with(params, {
    
    #update screening strategy
    sim_data <- assign_screening_strategy(params, sim_data)
    
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
    
    #update test trackers
    sim_data <- sim_data %>% mutate(
      num_colonoscopies = case_when(colonoscopy == 1 ~ as.integer(num_colonoscopies + 1), colonoscopy == 0 ~ num_colonoscopies),
      num_sigmoidoscopies = case_when(sigmoidoscopy == 1 ~ as.integer(num_sigmoidoscopies + 1), sigmoidoscopy == 0 ~ num_sigmoidoscopies),
      num_FIT = case_when(FIT == 1 ~ as.integer(num_FIT + 1), FIT == 0 ~ num_FIT),
      num_FOBT = case_when(FOBT == 1 ~ as.integer(num_FOBT + 1), FOBT == 0 ~ num_FOBT)
    )
    
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


assign_screening_strategy <- function(params, sim_data){
  
  with(params, {
    
    #update individual screening program type
    out <- sim_data %>%
      mutate(individual_screening_type = case_when( (age-starting_age+1) >= start_program_year ~ as.character(screening_type),
                                                    TRUE ~ "none"))
    
    #update an individual's next possible screening age based on results of last colonoscopy
    out <- out %>%
      mutate(
        next_possible_colonoscopy_age = case_when(
          clinical_distant_CRC_trt == 1 ~ as.integer(0),
          negative_colonoscopy == 1 ~ as.integer(next_possible_colonoscopy_age + 10),
          low_risk_polyp_identified == 1 ~ as.integer(next_possible_colonoscopy_age + 5),
          high_risk_polyp_identified == 1 ~ as.integer(next_possible_colonoscopy_age + 3),
          clinical_localized_CRC_trt == 1 | clinical_regional_CRC_trt == 1 ~ as.integer(age + 4),
          clinical_localized_CRC_trt == 0 & clinical_regional_CRC_trt == 0 &
            negative_colonoscopy == 0 & low_risk_polyp_identified == 0 & high_risk_polyp_identified == 0 &
            next_possible_colonoscopy_age == age - 1 ~ as.integer(next_possible_colonoscopy_age + 1),
          TRUE ~ next_possible_colonoscopy_age
        ),
        next_possible_sigmoidoscopy_age = case_when(
          !(individual_screening_type %in% c('FOBT + flex sig', 'flex sig')) ~ as.integer(0),
          hx_of_polyps == 1 | hx_of_CRC == 1 ~ as.integer(0),
          sigmoidoscopy == 1 & hx_of_polyps == 0 & hx_of_CRC == 0 &
            next_possible_sigmoidoscopy_age == age - 1 ~ as.integer(next_possible_sigmoidoscopy_age + 5),
          sigmoidoscopy == 0 & hx_of_polyps == 0 & hx_of_CRC == 0 &
            next_possible_sigmoidoscopy_age == age - 1  ~ as.integer(next_possible_sigmoidoscopy_age + 1),
          TRUE ~ next_possible_sigmoidoscopy_age
        ),
        hx_of_polyps = case_when(
          hx_of_polyps == 1 | low_risk_polyp_identified == 1 | high_risk_polyp_identified == 1 ~ as.integer(1),
          TRUE ~ as.integer(0)
        ),
        hx_of_CRC = case_when(
          hx_of_CRC == 1 | clinical_localized_CRC == 1 | clinical_regional_CRC == 1 | clinical_distant_CRC == 1 ~ as.integer(1),
          TRUE ~ as.integer(0)
        )
      )
    
    return(out)
    
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
        
      ) %>% dplyr::select(-risk)
    
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
          dead == 0 & (clinical_localized_CRC == 1 | clinical_regional_CRC == 1) ~ QALYs_gained + (QALY_local_regional_cancer*((1-discount_rate)^(age-starting_age))),
          dead == 0 & clinical_distant_CRC == 1 ~ QALYs_gained + (QALY_distant_cancer*((1-discount_rate)^(age-starting_age))),
          TRUE ~ QALYs_gained + (1*((1-discount_rate)^(age-starting_age)))
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
    ) %>% dplyr::select(-temp_dead)
    
    return(out)
    
  })
  
}


apply_screening <- function(params, sim_data){
  
  with(params, {

    out <- sim_data %>% mutate(
      
      #evaluate whether patient undergoes FOBT
      FOBT = case_when(age != next_possible_colonoscopy_age | hx_of_polyps == 1 | hx_of_CRC == 1 ~ as.integer(0),
                       individual_screening_type %in% c('FOBT', 'FOBT + flex sig') & dead == 0 & age <= max_screening_age &
                         (normal_mucosa == 1 | low_risk_polyp == 1 | high_risk_polyp == 1 | 
                            preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 | preclinical_distant_CRC == 1) ~ 
                         rbinom(n(), size=1, prob=adherence_FOBT),
                       TRUE ~ as.integer(0)
                       ),
      #evaluate whether patient undergoes FIT
      FIT = case_when(age != next_possible_colonoscopy_age | hx_of_polyps == 1 | hx_of_CRC == 1 ~ as.integer(0),
                      individual_screening_type == 'FIT' & dead == 0 & age <= max_screening_age &
                        (normal_mucosa == 1 | low_risk_polyp == 1 | high_risk_polyp == 1 | 
                          preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 | preclinical_distant_CRC == 1) ~ 
                        rbinom(n(), size=1, prob=adherence_FIT),
                      TRUE ~ as.integer(0)
      ),
      #evaluate whether patient undergoes flex sig
      sigmoidoscopy = case_when(age != next_possible_colonoscopy_age | hx_of_polyps == 1 | hx_of_CRC == 1 ~ as.integer(0),
                                individual_screening_type %in% c('FOBT + flex sig', 'flex sig') & age == next_possible_sigmoidoscopy_age & 
                                  dead == 0 & age <= max_screening_age &
                                  (normal_mucosa == 1 | low_risk_polyp == 1 | high_risk_polyp == 1 | 
                                     preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 | preclinical_distant_CRC == 1) ~
                                  rbinom(n(), size=1, prob=adherence_sigmoidoscopy),
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
      FIT_positive = case_when(FIT == 1 & normal_mucosa == 1 ~ rbinom(n(), size=1, prob=1-FIT_spec),
                               FIT == 1 & low_risk_polyp == 1 ~ rbinom(n(), size=1, prob=FIT_sens_low_risk_polyp),
                               FIT == 1 & high_risk_polyp == 1 ~ rbinom(n(), size=1, prob=FIT_sens_high_risk_polyp),
                               FIT == 1 & (preclinical_early_CRC == 1 | preclinical_regional_CRC == 1 | preclinical_distant_CRC == 1) ~ 
                                 rbinom(n(), size=1, prob=FIT_sens_CRC),
                               TRUE ~ as.integer(0)),
      sigmoidoscopy_positive = case_when(sigmoidoscopy == 1 & normal_mucosa == 1 ~ rbinom(n(), size=1, prob=1-sigmoidoscopy_spec),
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
      colonoscopy = case_when(age == next_possible_colonoscopy_age & 
                                (individual_screening_type == 'colonoscopy' | hx_of_polyps == 1 | hx_of_CRC == 1) & 
                                dead == 0 & age <= max_screening_age &
                                clinical_localized_CRC == 0 & clinical_regional_CRC == 0 & clinical_distant_CRC == 0 ~ 
                                rbinom(n(), size=1, prob=adherence_colonoscopy),
                              
                              #Other strategies
                              FOBT_positive == 1 | sigmoidoscopy_positive == 1 | FIT_positive == 1 ~ 
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
    
    out <- out %>% mutate(
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
      
      #indicate negative colonoscopy
      negative_colonoscopy = case_when(
        colonoscopy == 1 & dead == 0 & normal_mucosa == 1 ~ as.integer(1),
        TRUE ~ as.integer(0)
      ),
      
      #indicate low risk polyp identified
      low_risk_polyp_identified = case_when(
        colonoscopy == 1 & dead == 0 & low_risk_polyp == 1 ~ rbinom(n(), size=1, prob=colonoscopy_sigmoidoscopy_sens_low_risk_polyp),
        TRUE ~ as.integer(0)
      ),
      
      #indicate high risk polyp identified
      high_risk_polyp_identified = case_when(
        colonoscopy == 1 & dead == 0 & high_risk_polyp == 1 ~ rbinom(n(), size=1, prob=colonoscopy_sigmoidoscopy_sens_high_risk_polyp),
        TRUE ~ as.integer(0)
      ),
      
      #remove polyps when undergoing colonoscopy
      normal_mucosa = case_when(
        low_risk_polyp_identified == 1 ~ as.integer(1),
        high_risk_polyp_identified == 1 ~ as.integer(1),
        TRUE ~ normal_mucosa
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
    
    #update negative colonoscopy indicator
    out <- out %>% mutate(
      
      negative_colonoscopy = case_when(
        negative_colonoscopy == 1 ~ as.integer(1),
        colonoscopy == 1 & ((low_risk_polyp_identified == 0 & low_risk_polyp == 1) | (high_risk_polyp_identified == 0 & high_risk_polyp == 1) |
          (clinical_localized_CRC == 0 & preclinical_early_CRC == 1) | (clinical_regional_CRC == 0 & preclinical_regional_CRC == 1) |
          (clinical_distant_CRC == 0 & preclinical_distant_CRC == 1)) ~ as.integer(1),
        TRUE ~ as.integer(0)
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
      normal_mucosa_to_low_risk_polyp_prob = predict(normal_mucosa_to_low_risk_polyp_model, newdata = data.frame(age=age)),
      low_risk_polyp = case_when(
        dead == 0 & low_risk_polyp == 1 & high_risk_polyp == 1 ~ as.integer(0),
        dead == 0 & normal_mucosa == 1 ~ rbinom(n(), size=1, prob=normal_mucosa_to_low_risk_polyp_prob),
        TRUE ~ low_risk_polyp
      ),
      normal_mucosa = case_when(
        dead == 0 & normal_mucosa == 1 & low_risk_polyp == 1 ~ as.integer(0),
        TRUE ~ normal_mucosa
      )
    ) %>% dplyr::select(-normal_mucosa_to_low_risk_polyp_prob)
    
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
      cost_attained = cost_attained + ((colonoscopy_cost*colonoscopy + FOBT_cost*FOBT + FIT_cost*FIT + 
        sigmoidoscopy_cost*sigmoidoscopy + colonoscopy_perforation_tx_cost*perforation + 
        local_cancer_tx_cost*clinical_localized_CRC_trt + regional_cancer_tx_cost*clinical_regional_CRC_trt +
        distant_cancer_tx_cost*clinical_distant_CRC_trt + surveillance_cost*surveillance)*((1-discount_rate)^(age-starting_age)))
    )
    
    return(out)
    
  })
  
}