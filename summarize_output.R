#load dependencies
source('./src/dependencies.R')

# load files

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

none <- loadRData("./output/none_1e+07_50.Rdata")
FOBT <- loadRData("./output/FOBT_1e+07_50.Rdata")
FOBT_flex_sig <- loadRData("./output/FOBT + flex sig_1e+07_50.Rdata")
flex_sig <- loadRData("./output/flex sig_1e+07_50.Rdata")
FIT <- loadRData("./output/FIT_1e+07_50.Rdata")
colonoscopy <- loadRData("./output/colonoscopy_1e+07_50.Rdata")

none <- loadRData("./output/none_1e+07_5.Rdata")
FOBT <- loadRData("./output/FOBT_1e+07_5.Rdata")
FOBT_flex_sig <- loadRData("./output/FOBT + flex sig_1e+07_5.Rdata")
flex_sig <- loadRData("./output/flex sig_1e+07_5.Rdata")
FIT <- loadRData("./output/FIT_1e+07_5.Rdata")
colonoscopy <- loadRData("./output/colonoscopy_1e+07_5.Rdata")

none <- loadRData("./output/none_1e+07_1.Rdata")
FOBT <- loadRData("./output/FOBT_1e+07_1.Rdata")
FOBT_flex_sig <- loadRData("./output/FOBT + flex sig_1e+07_1.Rdata")
flex_sig <- loadRData("./output/flex sig_1e+07_1.Rdata")
FIT <- loadRData("./output/FIT_1e+07_1.Rdata")
colonoscopy <- loadRData("./output/colonoscopy_1e+07_1.Rdata")

# summaries

summaries <- 
  data.frame(
    
    'none' =
      c(none$num_colonoscopies %>% mean,
        none$num_sigmoidoscopies %>% mean,
        none$num_FOBT %>% mean,
        none$num_FIT %>% mean,
        none$cost_attained %>% mean,
        none$QALYs_gained %>% mean,
        none$dead %>% {sum(.)/length(.)},
        none$CRC_death %>% {sum(.)/length(.)}),
    
    'FOBT' =
      c(FOBT$num_colonoscopies %>% mean,
        FOBT$num_sigmoidoscopies %>% mean,
        FOBT$num_FOBT %>% mean,
        FOBT$num_FIT %>% mean,
        FOBT$cost_attained %>% mean,
        FOBT$QALYs_gained %>% mean,
        FOBT$dead %>% {sum(.)/length(.)},
        FOBT$CRC_death %>% {sum(.)/length(.)}),
    
    'FOBT + flex sig' = 
      c(FOBT_flex_sig$num_colonoscopies %>% mean,
        FOBT_flex_sig$num_sigmoidoscopies %>% mean,
        FOBT_flex_sig$num_FOBT %>% mean,
        FOBT_flex_sig$num_FIT %>% mean,
        FOBT_flex_sig$cost_attained %>% mean,
        FOBT_flex_sig$QALYs_gained %>% mean,
        FOBT_flex_sig$dead %>% {sum(.)/length(.)},
        FOBT_flex_sig$CRC_death %>% {sum(.)/length(.)}),
    
    'flex sig' = 
      c(flex_sig$num_colonoscopies %>% mean,
        flex_sig$num_sigmoidoscopies %>% mean,
        flex_sig$num_FOBT %>% mean,
        flex_sig$num_FIT %>% mean,
        flex_sig$cost_attained %>% mean,
        flex_sig$QALYs_gained %>% mean,
        flex_sig$dead %>% {sum(.)/length(.)},
        flex_sig$CRC_death %>% {sum(.)/length(.)}),
    
    'FIT' =
      c(FIT$num_colonoscopies %>% mean,
        FIT$num_sigmoidoscopies %>% mean,
        FIT$num_FOBT %>% mean,
        FIT$num_FIT %>% mean,
        FIT$cost_attained %>% mean,
        FIT$QALYs_gained %>% mean,
        FIT$dead %>% {sum(.)/length(.)},
        FIT$CRC_death %>% {sum(.)/length(.)}),
    
    'colonoscopy' =
      c(colonoscopy$num_colonoscopies %>% mean,
        colonoscopy$num_sigmoidoscopies %>% mean,
        colonoscopy$num_FOBT %>% mean,
        colonoscopy$num_FIT %>% mean,
        colonoscopy$cost_attained %>% mean,
        colonoscopy$QALYs_gained %>% mean,
        colonoscopy$dead %>% {sum(.)/length(.)},
        colonoscopy$CRC_death %>% {sum(.)/length(.)})
  )

rownames(summaries) <- c('av_num_colonoscopies', 'av_num_sigmoidoscopies', 'av_num_FOBT', 'av_num_FIT',
                         'cost attained', 'QALYs gained', 'proportion dead', 'proportion with CRC death')

summaries


# CRC mortality decrease compared to no screening

( (none$CRC_death %>% {sum(.)/length(.)})-(none$CRC_death %>% {sum(.)/length(.)}) ) / (none$CRC_death %>% {sum(.)/length(.)})
( (none$CRC_death %>% {sum(.)/length(.)})-(FOBT$CRC_death %>% {sum(.)/length(.)}) ) / (none$CRC_death %>% {sum(.)/length(.)})
( (none$CRC_death %>% {sum(.)/length(.)})-(FOBT_flex_sig$CRC_death %>% {sum(.)/length(.)}) ) / (none$CRC_death %>% {sum(.)/length(.)})
( (none$CRC_death %>% {sum(.)/length(.)})-(flex_sig$CRC_death %>% {sum(.)/length(.)}) ) / (none$CRC_death %>% {sum(.)/length(.)})
( (none$CRC_death %>% {sum(.)/length(.)})-(FIT$CRC_death %>% {sum(.)/length(.)}) ) / (none$CRC_death %>% {sum(.)/length(.)})
( (none$CRC_death %>% {sum(.)/length(.)})-(colonoscopy$CRC_death %>% {sum(.)/length(.)}) ) / (none$CRC_death %>% {sum(.)/length(.)})


# create ICER table

none_cost <- none$cost_attained %>% mean()
FOBT_cost <- FOBT$cost_attained %>% mean()
FOBT_flex_sig_cost <- FOBT_flex_sig$cost_attained %>% mean()
flex_sig_cost <- flex_sig$cost_attained %>% mean()
FIT_cost <- FIT$cost_attained %>% mean()
colonoscopy_cost <- colonoscopy$cost_attained %>% mean()

none_effect <- none$QALYs_gained %>% mean()
FOBT_effect <- FOBT$QALYs_gained %>% mean()
FOBT_flex_sig_effect <- FOBT_flex_sig$QALYs_gained %>% mean()
flex_sig_effect <- flex_sig$QALYs_gained %>% mean()
FIT_effect <- FIT$QALYs_gained %>% mean()
colonoscopy_effect <- colonoscopy$QALYs_gained %>% mean()

CE_df <- data.frame(strategy = c('none', 'FOBT', 'FOBT + flex sig', 'flex sig', 'FIT', 'colonoscopy'),
                    cost = c(none_cost, FOBT_cost, FOBT_flex_sig_cost, flex_sig_cost, FIT_cost, colonoscopy_cost),
                    effect = c(none_effect, FOBT_effect, FOBT_flex_sig_effect, flex_sig_effect, FIT_effect, colonoscopy_effect))

strongly_dominated <- sapply(1:nrow(CE_df), function(x) {
  any(CE_df$cost[x] > CE_df$cost & CE_df$effect[x] < CE_df$effect)
})

CE_df <- CE_df[!strongly_dominated,] %>% 
  {.[order(.$cost),]} %>%
  mutate(
    incremental_cost = cost-lag(cost),
    incremental_effect = effect-lag(effect),
    ICER = incremental_cost/incremental_effect
  )

while (is.unsorted(CE_df$ICER, na.rm=T)) {
  
  CE_df <- CE_df %>% 
    mutate(weakly_dominated = case_when(ICER > lead(ICER) ~ 1, TRUE ~ 0)) %>%
    filter(row_number() != which(weakly_dominated == 1)[1]) %>%
    select(-weakly_dominated) %>%
    mutate(
      incremental_cost = cost-lag(cost),
      incremental_effect = effect-lag(effect),
      ICER = incremental_cost/incremental_effect
    )
  
}

CE_df


# plot frontier

CE_df_plot <- data.frame(strategy = c('none', 'FOBT', 'FOBT + flex sig', 'flex sig', 'FIT', 'colonoscopy'),
                         cost = c(none_cost, FOBT_cost, FOBT_flex_sig_cost, flex_sig_cost, FIT_cost, colonoscopy_cost),
                         effect = c(none_effect, FOBT_effect, FOBT_flex_sig_effect, flex_sig_effect, FIT_effect, colonoscopy_effect))

CE_df_plot %>% ggplot(data=., aes(x=effect, y=cost, col=strategy)) + geom_point()
