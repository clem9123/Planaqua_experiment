#################### Model 14
# 2016 avec Lake variable et pas t0

inits <- function(){
  list(t0 = runif(1, 2010, 2016), 
       K = runif(4,0,10), 
       Kla = runif(16,0,10),
       sdKla = runif(4,0,20),
       Linf = runif(4,130,250), 
       Linfla = runif(16,0,10),
       sdLla = runif(4,0,100),
       sd = runif(1,0,20))}

parameters = c("K","Linf","sd","Kla","Linfla", "sdKla", "sdLla","t0")

model_g_Ktrla_Ltrla <- function ()
{
  #Priors
  sd ~ dunif(0,20)
  t0 ~ dunif(2010,2016)
  for(tr in 1:4){
    Linf[tr] ~ dunif (130,250)
    K[tr] ~ dunif(0,10)
    sdKla[tr] ~ dunif(0,20)
    sdLla[tr] ~ dunif(0,100)
  }
  for (la in 1:16){
    Linfla[la] ~ dnorm (Linf[Lt[la,2]], 1/(sdLla[Lt[la,2]]^2))
    Kla[la] ~ dnorm (K[Lt[la,2]], 1/(sdKla[Lt[la,2]]^2))
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linfla[Lake[i]]*(1-exp(-Kla[Lake[i]]*(2015+t-t0))), 1/(sd^2))
    } # t
  } # i
} # func

Model_g_Ktrla_Ltrla_2016 <- jags.parallel(data = jags.data,
                                          inits = inits,
                                          parameters.to.save = parameters,
                                          model.file = model_g_Ktrla_Ltrla,
                                          n.chains = 4,
                                          n.iter = 5000)

#################### Model 14
# 2016 avec Lake variable et pas t0

inits <- function(){
  list(t0 = runif(1, 2010, 2016), 
       K = runif(4,0,10), 
       Kla = runif(16,0,10),
       sdKla = runif(1,0,20),
       Linf = runif(4,130,250), 
       Linfla = runif(16,0,10),
       sdLla = runif(1,0,100),
       sd = runif(1,0,20))}

parameters = c("K","Linf","sd","Kla","Linfla", "sdKla", "sdLla","t0")

model_g_Ktrla_Ltrla <- function ()
{
  #Priors
  sd ~ dunif(0,20)
  t0 ~ dunif(2010,2016)
  sdKla ~ dunif(0,20)
  sdLla ~ dunif(0,100)
  for(tr in 1:4){
    Linf[tr] ~ dunif (130,250)
    K[tr] ~ dunif(0,10)
    
  }
  for (la in 1:16){
    Linfla[la] ~ dnorm (Linf[Lt[la,2]], 1/(sdLla^2))
    Kla[la] ~ dnorm (K[Lt[la,2]], 1/(sdKla^2))
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linfla[Lake[i]]*(1-exp(-Kla[Lake[i]]*(2015+t-t0))), 1/(sd^2))
    } # t
  } # i
} # func

Model_g_Ktrla_Ltrla_2016 <- jags.parallel(data = jags.data,
                                          inits = inits,
                                          parameters.to.save = parameters,
                                          model.file = model_g_Ktrla_Ltrla,
                                          n.chains = 4,
                                          n.iter = ni)

