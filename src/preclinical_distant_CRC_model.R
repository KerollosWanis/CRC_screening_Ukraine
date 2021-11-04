data.frame(
  age=c(50,55,60,65,70,75), 
  probability=c(0.0004,
                0.0002,
                0.0004,
                0.0005,
                0.0006,
                0.0008)
) %>% lm(data=., probability ~ age + I(age^2) + I(log(age))) 