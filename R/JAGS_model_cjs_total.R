################# data


tr <- BDD_a %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
mooved <- tr[duplicated(tr$Tag_id),]$Tag_id
tr <- tr %>% filter ( !duplicated(tr$Tag_id))
tr <- tr %>% mutate (Treatment = ifelse(Tag_id %in% mooved, NA, Treatment))

la <- BDD_a %>% group_by(Tag_id, Lake) %>% summarize() %>% ungroup()
mooved <- la[duplicated(la$Tag_id),]$Tag_id
la <- la %>% filter ( !duplicated(la$Tag_id))
la <- la %>% mutate (Lake = ifelse(Tag_id %in% mooved, NA, Lake))

t <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from =Lake)
t <- merge(t,tr)
t <- merge(t,la)

t <- data.frame(t$Tag_id, apply(t[2:8],2, function(x) ifelse(is.na(x),0,1)), t$Treatment, t$Lake)
t <- na.omit(t)

set.seed(10)
t <- t[sample(1:nrow(t), 100),]

s <- s %>% filter (Tag_id %in% t$t.Tag_id)
t <- t %>% arrange(t.Tag_id)
s <- s %>% arrange(Tag_id)

CH <- as.matrix(t[2:7])

get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

z.inits <- function(ch){
  state <- ch
  state[state==0] <- 1
  get.first <- function(x) min(which(x!=0))
  f <- apply(ch, 1, get.first)
  for (i in 1:nrow(ch)){
    state[i,1:(f[i])] <- NA
  }
  return(state)
}

jags.data <- list(y = CH,
                  size = as.matrix(s[2:7]),
                  Treatment = t$t.Treatment,
                  f = apply(CH, 1, get.first), 
                  nind = nrow(CH), 
                  noccas = ncol(CH),
                  ni = 10000,
                  zi = z.inits(CH))

############################ Model 1 phi(~1)p(~1)

inits <- function(){
  list(phi = runif(0,1), p = runif(0,1), z = zi)}


parameters <- c("phi", "p")

cjs_0 <- function(){
  # Likelihood
  for (i in 1:nind){
    # State at first capture is alive
    z[i,f[i]] <- 1 
    for (t in (f[i]+1):noccas){
      # determine the state alive/dead
      z[i,t] ~ dbern(phi * z[i,t-1]) 
      # determine the capture status
      y[i,t] ~ dbern(p * z[i,t]) 
    }
  }
  # Priors
  phi ~ dunif(0,1)
  p ~ dunif(0,1)
}

CJS_0 <- jags.parallel(data = jags.data,
                       inits = inits,
                       parameters.to.save = parameters,
                       model.file = cjs_0,
                       n.chains = 4,
                       n.iter = ni)

save(CJS_0, file = "object/CJS_0.RData")

##################################### Model 2 phi(~t+i)p(~t+i)

inits <- function() {
  list(phi = 0.2, Phi = matrix(0.5, ncol = noccas-1, nrow = nind) , p = 0.1, P = matrix(0.5, ncol = noccas-1, nrow = nind), z = zi)}

parameters = c("phi", "p")

cjs_1 <- function()
{
  # Priors
  for (i in 1:nind){
    for (t in 1:(noccas-1)){
      Phi[i,t] ~ dnorm(phi,0.1) # Prior individual and time dependant survival
      P[i,t] ~ dnorm (p,0.1) # Prior individual and time dependant capture
    } #t
  } #i        
  p ~ dunif(0, 1) # Prior for mean recapture
  phi ~ dunif(0, 1) 
  
  # Likelihood 
  for (i in 1:nind){
    
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      # survival process
      
      
      # State process
      z[i,t] ~ dbern(Phi[i,t-1] * z[i,t-1])
      
      # Observation process
      y[i,t] ~ dbern(P[i,t-1] * z[i,t])
    } #t
  } #i
}

CJS_1 <- jags.parallel(data = jags.data, 
              inits = inits, 
              parameters.to.save = parameters,
              model.file = cjs_1,
              n.chains = 4,
              n.iter = ni)

save(CJS_1, file = "object/CJS_1.RData")

################# Model 1 phi(~Treatment)p(~1)

inits <- function(){
  list(p = runif(1,0,1), phi = runif(4,0,1), z = zi)
}

parameters = c("phi","p")


cjs_treatment <- function(){
  # likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(phi[Treatment[i]]*z[i, t-1])
      y[i,t] ~ dbern(p*z[i,t])
    }
  }
  # Priors
  p~ dunif(0,1)
  for (tr in 1:4){
    phi[tr] ~ dunif(0,1)
  } 
}

CJS_treatment <- jags.parallel(data = jags.data,
                               inits = inits,
                               parameters.to.save = parameters,
                               model.file = cjs_treatment,
                               n.chains = 4,
                               n.iter = ni)
#print(CJS_treatment)
#autocorr.plot(CJS_treatment, ask = F)
#traceplot(CJS_treatment, ask = F)

save(CJS_treatment, file = "R/object/CJS_treatment.RData")

############### Model 2 phi(~Treatment + age)p(~age)

inits <- function(){
  list(p = runif(6,0,1), phi = matrix(ncol = 6, runif(24,0,1)), z = zi)
}

parameters = c("phi","p")

cjs_treatment_age <- function() {
  #likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(phi[Treatment[i], t-f[i]]*z[i,t-1])
      y[i,t] ~ dbern(p[t-f[i]]*z[i,t])
    }
  }
  #prior
  for (age in 1:6){
    p[age] ~ dunif(0,1)
    for (tr in 1:4){
      phi[tr, age] ~ dunif(0,1)
    }
  }
}

CJS_treatment_age <- jags.parallel(data = jags.data,
                                   inits = inits,
                                   parameters.to.save = parameters,
                                   model.file = cjs_treatment_age,
                                   n.chains = 4,
                                   n.iter = ni)
#print(CJS_treatment_age)
#autocorr.plot(CJS_treatment_age, ask = F)
#traceplot(CJS_treatment_age, ask = F)

save(CJS_treatment_age, file = "R/object/CJS_treatment_age.RData")

################## Model 3 phi(~Treatment + maturity)p(~ maturity)

inits <- function(){
  list(pm = runif(2,0,1), phim = matrix(ncol = 2, runif(8,0,1)), z = zi)
}

parameters = c("phim","pm")

cjs_treatment_maturity <- function() {
  #likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(phi[Treatment[i], t-f[i]]*z[i,t-1])
      y[i,t] ~ dbern(p[t-f[i]]*z[i,t])
    }
  }
  #prior
  for (m in 1:2){
    pm[m] ~ dunif(0,1)
    for (tr in 1:4){
      phim[tr, m] ~ dunif(0,1)
    }
  }
  #constraint
  p[1] <- pm[1]
  for (tr in 1:4){
    phi[tr, 1] <- phim[tr, 1]
  }
  for (age in 2:6){
    p[age] <- pm[2]
    for (tr in 1:4){
      phi[tr, age] <- phim[tr, 2]
    }
  }
}

CJS_treatment_maturity <- jags.parallel(data = jags.data,
                                        inits = inits,
                                        parameters.to.save = parameters,
                                        model.file = cjs_treatment_maturity,
                                        n.chains = 4,
                                        n.iter = ni)
#print(CJS_treatment_maturity)
#autocorr.plot(CJS_treatment_maturity, ask = F)
#traceplot(CJS_treatment_maturity, ask = F)

save(CJS_treatment_maturity, file = "R/object/CJS_treatment_maturity.RData")

