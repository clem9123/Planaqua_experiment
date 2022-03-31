library(ggplot2)
library(R2jags)

n <- 1000
Height <- rnorm(n, mean = 170, sd = 20)

ggplot(data.frame(Height))+
  geom_histogram(aes(x = Height))

Height.jags = list(Height = Height, n = length(Height))
inits = list(list(mean = 170, sd = 30), list(mean = 100, sd = 50))
parameters = c("mean", "sd")

cat(file = "jags_models/model_test.txt", "
model
{
  # Likelihood function
  for (i in 1 : n){
    Height[i] ~ dnorm(mean, sd^2)
    }
# Priors
mean ~ dnorm(0,0.0001)
sd ~ dnorm(0,0.0001)
}
")

model_test <- function()
{
  # Likelihood function
  for (i in 1 : n){
    Height[i] ~ dnorm(mean, sd^2)
    }
# Priors
mean ~ dnorm(0,0.0001)
sd ~ dnorm(0,0.0001)
}



test <- jags(data = Height.jags, 
                inits = inits, 
                parameters.to.save = parameters, 
                model.file = model_test, 
                n.chains = 2,
                n.iter = 2000, 
                n.burnin = 1000)

print(test, digits = 3)
traceplot(test)


# visualisation de mes résulats d'estimate

est <- c(as.mcmc(test))
mean <- c(data.frame(est[1])$mean, data.frame(est[2])$mean)

ggplot(data.frame(mean))+
  geom_density(aes(x = mean))

sd <- c(data.frame(est[1])$sd, data.frame(est[2])$sd)

ggplot(data.frame(sd))+
  geom_density(aes(x = sd))

predict = rnorm(1000,mean(mean), 1/mean(sd))
ggplot(data.frame(predict))+
  geom_density(aes(x = predict))

ggplot(data.frame(Height))+
  geom_density(aes(x = Height))+
  geom_density(aes(x = predict), color = "red")

t = rnorm(1000,0,1/0.001)
ggplot(data.frame(t))+
  geom_density(aes(x = t))

# C'est bon je pense que j'ai compris le principe de Jags :D
# Maintenant il va falloir comprendre les priors et la rédaction de model