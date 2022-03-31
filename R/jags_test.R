library(ggplot2)
library(R2jags)

n <- 1000
Height <- rnorm(n , mean = 170, sd = 10)

ggplot(data.frame(Height))+
  geom_histogram(aes(x = Height))

Height.jags = list(Height = Height, n = length(Height))
inits = list(list(mean = 170, sd = 12), list(mean = 100, sd = 1))
parameters = c("mean", "sd")

cat(file = "jags_models/model_test.txt", "
model
{
  # Likelihood function
  for (i in 1 : n){
    Height[i] ~ dnorm(mean,sd)
    }
# Priors
mean ~ dnorm(0,1)
sd ~ dnorm(1,0.01)
}
")


test <- jags(data = Height.jags, 
                inits = inits, 
                parameters.to.save = parameters, 
                model.file = "jags_models/model_test.txt", 
                n.chains = 2,
                n.iter = 2000, 
                n.burnin = 1000)

print(test, digits = 3)
traceplot(test)
mean <- test$sims.list$mean

ggplot(data.frame(mean))+
  geom_density(aes(x = mean))
