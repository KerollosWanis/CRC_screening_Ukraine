# load files

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
none <- loadRData("./output/none_1e+07_50.Rdata")
FOBT <- loadRData("./output/FOBT_1e+07_50.Rdata")
FOBT_flex_sig <- loadRData("./output/FOBT + flex sig_1e+07_50.Rdata")
colonoscopy <- loadRData("./output/colonoscopy_1e+07_50.Rdata")


# summaries

c(none$cost_attained %>% mean,
  none$QALYs_gained %>% mean,
  none$colonoscopy %>% sum,
  none$dead %>% {sum(.)/length(.)},
  none$CRC_death %>% {sum(.)/length(.)})

c(FOBT$cost_attained %>% mean,
  FOBT$QALYs_gained %>% mean,
  FOBT$colonoscopy %>% sum,
  FOBT$dead %>% {sum(.)/length(.)},
  FOBT$CRC_death %>% {sum(.)/length(.)})

c(FOBT_flex_sig$cost_attained %>% mean,
  FOBT_flex_sig$QALYs_gained %>% mean,
  FOBT_flex_sig$colonoscopy %>% sum,
  FOBT_flex_sig$dead %>% {sum(.)/length(.)},
  FOBT_flex_sig$CRC_death %>% {sum(.)/length(.)})

c(colonoscopy$cost_attained %>% mean,
  colonoscopy$QALYs_gained %>% mean,
  colonoscopy$colonoscopy %>% sum,
  colonoscopy$dead %>% {sum(.)/length(.)},
  colonoscopy$CRC_death %>% {sum(.)/length(.)})


# CRC mortality decrease compared to no screening

( (none$CRC_death %>% {sum(.)/length(.)})-(none$CRC_death %>% {sum(.)/length(.)}) ) / (none$CRC_death %>% {sum(.)/length(.)})
( (none$CRC_death %>% {sum(.)/length(.)})-(FOBT$CRC_death %>% {sum(.)/length(.)}) ) / (none$CRC_death %>% {sum(.)/length(.)})
( (none$CRC_death %>% {sum(.)/length(.)})-(FOBT_flex_sig$CRC_death %>% {sum(.)/length(.)}) ) / (none$CRC_death %>% {sum(.)/length(.)})
( (none$CRC_death %>% {sum(.)/length(.)})-(colonoscopy$CRC_death %>% {sum(.)/length(.)}) ) / (none$CRC_death %>% {sum(.)/length(.)})


# create ICER table

none_cost <- none$cost_attained %>% mean()
FOBT_cost <- FOBT$cost_attained %>% mean()
FOBT_flex_sig_cost <- FOBT_flex_sig$cost_attained %>% mean()
colonoscopy_cost <- colonoscopy$cost_attained %>% mean()

none_effect <- none$QALYs_gained %>% mean()
FOBT_effect <- FOBT$QALYs_gained %>% mean()
FOBT_flex_sig_effect <- FOBT_flex_sig$QALYs_gained %>% mean()
colonoscopy_effect <- colonoscopy$QALYs_gained %>% mean()

CE_df <- data.frame(strategy = c('none', 'FOBT', 'FOBT + flex sig', 'colonoscopy'),
                    cost = c(none_cost, FOBT_cost, FOBT_flex_sig_cost, colonoscopy_cost),
                    effect = c(none_effect, FOBT_effect, FOBT_flex_sig_effect, colonoscopy_effect))

CE_df <- CE_df[order(CE_df$cost),] %>% mutate(
  incremental_cost = cost-lag(cost),
  incremental_effect = effect-lag(effect),
  ICER = incremental_cost/incremental_effect
) %>% filter(ICER > 0 | is.na(ICER)) %>% mutate(
  incremental_cost = cost-lag(cost),
  incremental_effect = effect-lag(effect),
  ICER = incremental_cost/incremental_effect
) %>% filter(ICER > 0 | is.na(ICER)) %>% mutate(
  incremental_cost = cost-lag(cost),
  incremental_effect = effect-lag(effect),
  ICER = incremental_cost/incremental_effect
)

CE_df


# plot frontier

CE_df_plot <- data.frame(strategy = c('none', 'FOBT', 'FOBT + flex sig', 'colonoscopy'),
                         cost = c(none_cost, FOBT_cost, FOBT_flex_sig_cost, colonoscopy_cost),
                         effect = c(none_effect, FOBT_effect, FOBT_flex_sig_effect, colonoscopy_effect))

CE_df_plot %>% ggplot(data=., aes(x=effect, y=cost, col=strategy)) + geom_point()