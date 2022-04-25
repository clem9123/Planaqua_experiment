##################### Model phi(~taille)p(~taille)


# logit(phi) <- / ~ c*taille +d
# logit(p) <- / ~ a*taille + b

inits <- function(){
  list(a = runif(1,-0.2,0.2), b = runif(1,-0.2,0.2), z = zi)
}

parameters = c("a","b")

cjs_taille1 <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      p[i,t_bis] <- 0
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      logit(p[i,t]) <- a*size[i,t]+2
      logit(phi[i,t]) <- b*size[i,t]+2
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  a ~ dunif(-0.2,0.2)
  b ~ dunif(-0.2,0.2)
}

CJS_taille1 <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = cjs_taille1,
                             n.chains = 4,
                             n.iter = ni)
#print(CJS_taille)
#autocorr.plot(CJS_taille, ask = F)
#traceplot(CJS_taille, ask = F)

save(CJS_taille, file = "R/object/CJS_taille.RData")

##################### Model phi(~taille)p(~taille)


# logit(phi) <- / ~ c*taille +d
# logit(p) <- / ~ a*taille + b

inits <- function(){
  list(a_p = runif(1,-0.2,0.2), a_phi = runif(1,-0.2,0.2),
       a2_p = runif(1,-0.2,0.2), a2_phi = runif(1,-0.2,0.2),
       b_p = runif(1,-4,1), b_phi = runif(1,-4,1),
       z = zi)
}

parameters = c("a_p","a_phi")

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
      logit(phi[i,t]) <- a2_phi*size[i,t]^2 + a_phi*size[i,t] + b_phi # + b_phi[Treatment[i]]
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  a_p ~ dunif(-0.2,0.2)
  a_phi ~ dunif(-0.2,0.2)
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

##################### Model phi(~taille)p(~taille)


# logit(phi) <- / ~ c*taille +d
# logit(p) <- / ~ a*taille + b

inits <- function(){
  list(b_p = runif(1,-3,1), b_phi = runif(1,-5,-2),
       z = zi)
}

parameters = c("b_p","b_phi")

cjs_taille2 <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      p[i,t_bis] <- 0
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      logit(p[i,t]) <- 0.08*size[i,t]^2 - 0.0004 *size[i,t] + b_p # + b_p[Treatment[i]]
      logit(phi[i,t]) <- 0.04*size[i,t] + b_phi # + b_phi[Treatment[i]]
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  b_p ~ dunif(-3,1)
  b_phi ~ dunif(-5,-2)
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
       a2_p = -4, #a2_phi = -0.00025,
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
      logit(p[i,t]) <- a2_p*10^-4*size[i,t]^2 + a_p*10^-2*size[i,t] + b_p # + b_p[Treatment[i]]
      logit(phi[i,t]) <- a_phi*10^-2*size[i,t] + b_phi # + b_phi[Treatment[i]]
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  a_p ~ dunif(5,11)
  a_phi ~ dunif(1,40)
  a2_p ~ dunif(-7,-2.5)
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


cjs_taille2 <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      #logit(p) <- 1/(1+exp(-(a2_p*10^-4*size[i,t]^2 + a_p*size[i,t] + b_p))) # + b_p[Treatment[i]]
      logit(phi[i,t]) <- a_phi*size[i,t] + b_phi # + b_phi[Treatment[i]]
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(z[i,t]/(1+exp(-(a2_p*10^-4*size[i,t]^2 + a_p*size[i,t] + b_p))))
    }
  }
  # Priors
  a_p ~ dunif(-0.2,0.2)
  a_phi ~ dunif(-0.2,0.2)
  a2_p ~ dunif(-20,20)
  #a2_phi ~ dunif(-0.2,0.2)
  b_p ~ dunif(-4,2)
  b_phi ~ dunif(-7,-1)
}