data.frame(
  age=c(50,55,60,65,70), 
  probability=c(0.00836,
                0.0099,
                0.01156,
                0.0133,
                0.01521)
) %>% lm(data=., probability ~ age + I(age^2) + I(log(age))) 