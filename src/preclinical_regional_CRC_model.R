data.frame(
  age=c(50,55,60,65,70,75), 
  probability=c(0.0012,
                0.0023,
                0.0040,
                0.0053,
                0.0062,
                0.0072)
) %>% lm(data=., probability ~ age + I(age^2) + I(log(age))) 