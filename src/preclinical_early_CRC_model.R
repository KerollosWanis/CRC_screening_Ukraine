data.frame(
  age=c(50,55,60,65,70,75), 
  probability=c(0.0024, 
                0.0067, 
                0.0092, 
                0.011, 
                0.013, 
                0.014)
) %>% lm(data=., probability ~ age + I(age^2) + I(log(age))) 