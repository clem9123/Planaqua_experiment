
##################### Model phi(~taille)p(~taille)


# logit(phi) <- / ~ c*taille +d
# logit(p) <- / ~ a*taille + b

inits <- function(){
  list(a_p = runif(1,-0.2,0.2), a_phi = runif(1,-0.2,0.2),
       a2_p = runif(1,-0.2,0.2), #a2_phi = runif(1,-0.2,0.2),
       b_p = runif(1,-4,2), b_phi = runif(1,-7,-1),
       z = zi)
}

inits <- function(){
  list(a_p = 8, a_phi = 4,
       a2_p = -4, #a2_phi = -3,
       b_p = -2, b_phi = -2.5,
       z = zi)
}

parameters = c("a2_p","a_p","a_phi","b_p","b_phi")

cjs_taille2 <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      p[i,t_bis] <- 0
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      logit(p[i,t]) <- a2_p*size[i,t]^2 + a_p*size[i,t] + b_p # + b_p[Treatment[i]]
      logit(phi[i,t]) <- a_phi*size[i,t] + b_phi # + b_phi[Treatment[i]]
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  a_p ~ dunif(5,11)
  a_phi ~ dunif(1,40)
  a2_p ~ dunif(-7,0)
  #a2_phi ~ dunif(-0.2,0.2)
  b_phi ~ dunif(-20,10)
  b_p ~ dunif(-4,-2)
}


CJS_taille2 <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = cjs_taille2,
                             n.chains = 4,
                             n.iter = ni)
#print(CJS_taille)
#autocorr.plot(CJS_taille, ask = F)
#traceplot(CJS_taille, ask = F)

save(CJS_taille2, file = "R/object/CJS_taille2.RData")


##################### Model phi(~taille + Treatment)p(~taille)


# logit(phi) <- / ~ c*taille +d
# logit(p) <- / ~ a*taille + b

inits <- function(){
  list(a_p = runif(1,7,10), a_phi = runif(1,8,15),
       a2_p = runif(1,-7,-1), #a2_phi = runif(1,-0.2,0.2),
       b_p = runif(1,-4,-2), b_phi = runif(4,-5,1),
       z = zi)
}

inits <- function(){
  list(a2_p = -4, a_p = 8, b_p = -2,
       a_phi = 12, b_phi = rep(-8,4), #a2_phi = -3.3, 
       z = zi)
}

parameters = c("a2_p","a_p","a_phi","b_p","b_phi")

cjs_taille2 <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      p[i,t_bis] <- 0
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      logit(p[i,t]) <- a2_p*size[i,t]^2 + a_p*size[i,t] + b_p # + b_p[Treatment[i]]
      logit(phi[i,t]) <- a_phi*size[i,t] + b_phi[Treatment[i]]
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  a_p ~ dunif(5,11)
  a_phi ~ dunif(8,15)
  a2_p ~ dunif(-7,0)
  #a2_phi ~ dunif(-7,0)
  for (tr in 1:4){
    b_phi[tr] ~ dunif(-5,1)
  }
  b_p ~ dunif(-4,-2)
}


CJS_taille2 <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = cjs_taille2,
                             n.chains = 4,
                             n.iter = ni)
#print(CJS_taille)
#autocorr.plot(CJS_taille, ask = F)
#traceplot(CJS_taille, ask = F)

save(CJS_taille2, file = "R/object/CJS_taille2.RData")

######################################## Model taille discrete


inits <- function(){
  list(phi= runif(6,0.4,0.8), p = runif(6,0.2,0.8))
}

parameters = c("phi", "p")

cjs_group_taille <- function() {
  #likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(phi[size_group[i,t]]*z[i,t-1])
      y[i,t] ~ dbern(p[size_group[i,t]]*z[i,t])
    }
  }
  #prior
  for (gs in 1:6){
    p[gs] ~ dunif(0,1)
    phi[gs] ~ dunif(0,1)
  }
}

CJS_group_taille <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = cjs_group_taille,
                             n.chains = 4,
                             jags.seed = 143,
                             n.iter = ni)
