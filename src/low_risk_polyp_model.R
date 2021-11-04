data.frame(
  age=c(50,60,70), 
  probability=c(0.2,0.4,0.5)
) %>% lm(data=., probability ~ age + I(age^2)) 