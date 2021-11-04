data.frame(
  age=c(50,60,70,80), 
  mortality=c(0.05,
              0.09,
              0.16,
              0.21)
) %>% lm(data=., mortality ~ age + I(age^2) + I(log(age))) 