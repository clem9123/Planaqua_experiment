################# data
library(R2jags)
library(arules)
s <- s2

tr <- BDD_a %>% filter(Year != "2022") %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
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

t <- data.frame(t$Tag_id, apply(t[2:7],2, function(x) ifelse(is.na(x),0,1)), t$Treatment, t$Lake)
t <- na.omit(t)

#set.seed(632)
#t <- t[sample(1:nrow(t), 500),]

s <- s %>% filter (Tag_id %in% t$t.Tag_id)
t <- t %>% arrange(t.Tag_id)
s <- s %>% arrange(Tag_id)

size_break <- c(0,120,140,180,300)
s_group <- data.frame(apply(s[2:7], 2, 
    function(x) discretize(x,method = "fixed",
                           breaks = size_break, 
                           labels= c(1:(length(size_break)-1)))))
s_group [] <- apply(s_group, 2, as.numeric)

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
                  size = as.matrix(s[2:7]*10^-2),
                  size_group = s_group,
                  Treatment = t$t.Treatment,
                  Lake = t$t.Lake,
                  f = apply(CH, 1, get.first), 
                  nind = nrow(CH), 
                  noccas = ncol(CH),
                  ni = 10000,
                  n_size = length(size_break)-1,
                  zi = z.inits(CH))

############################ Model 1 phi(~1)p(~1)

inits <- function(){
  list(phi = runif(1,0,1), p = runif(1,0,1), z = zi)}


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



######################################## Model taille discrete


inits <- function(){
  list(phi= runif(n_size,0,1), p = runif(n_size,0,1), z =zi) #
}

parameters = c("phi", "p")

cjs_group_taille <- function() {
  #likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(phi[size_group[i,t]]*z[i,t-1])#
      y[i,t] ~ dbern(p[size_group[i,t]]*z[i,t])#
    }
  }
  #prior
  for (gs in 1:n_size){
    p[gs] ~ dunif(0,1)
    phi[gs] ~ dunif(0,1)
  }
}


CJS_group_taille <- jags.parallel(data = jags.data,
                                  inits = inits,
                                  parameters.to.save = parameters,
                                  model.file = cjs_group_taille,
                                  n.chains = 4,
                                  n.iter = ni)

save(CJS_group_taille, file = "R/object/CJS_group_taille.RData")

######################################## Model taille discrete


inits <- function(){
  list(p= runif(n_size,0.1,0.9), 
       phi = matrix(ncol = 4, runif(4*n_size,0.1,0.9)),
       z = zi)
}

parameters = c("phi", "p")

cjs_group_taille_tr <- function() {
  #likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(phi[size_group[i,t-1], Treatment[i]]*z[i,t-1])
      y[i,t] ~ dbern(p[size_group[i,t]]*z[i,t])
    }
  }
  #prior
  for (gs in 1:n_size){
      p[gs] ~ dunif(0,1)
    
    for (tr in 1:4){
      phi[gs,tr] ~ dunif(0,1)
    }
  }
}

CJS_group_taille_tr <- jags.parallel(data = jags.data,
                                       inits = inits,
                                       parameters.to.save = parameters,
                                       model.file = cjs_group_taille_tr,
                                       n.chains = 4,
                                       n.iter = ni)

save(CJS_group_taille_tr, file = "R/object/CJS_group_taille_tr.RData")

######################################## Model taille discrete


inits <- function(){
  list(p= array(rep(runif(1,0.1,0.9), n_size*noccas*16), c(n_size,noccas,16)), 
       phi = matrix(ncol = 4, runif(4*n_size,0.1,0.9)),
       z = zi)
}

parameters = c("phi", "p")

cjs_group_taille_tr_p <- function() {
  #likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(phi[size_group[i,t-1], Treatment[i]]*z[i,t-1])
      y[i,t] ~ dbern(p[size_group[i,t],t, Lake[i]]*z[i,t])
    }
  }
  #prior
  for (gs in 1:n_size){
    for (t in 1:noccas){
      
      for (l in 1:16){
        p[gs,t,l] ~ dunif(0,1)
      }
    }
    
    for (tr in 1:4){
      phi[gs,tr] ~ dunif(0,1)
    }
  }
}

CJS_group_taille_tr_p <- jags.parallel(data = jags.data,
                                  inits = inits,
                                  parameters.to.save = parameters,
                                  model.file = cjs_group_taille_tr_p,
                                  n.chains = 4,
                                  n.iter = ni)

save(CJS_group_taille_tr_p, file = "R/object/CJS_group_taille_tr_p.RData")

######################################## Model taille discrete


inits <- function(){
  list(p= array(rep(runif(1,0.1,0.9), n_size*(noccas-1)*16), c(n_size,(noccas-1),16)), 
       phi = array(rep(runif(1,0.1,0.9), n_size*(noccas-1)*4), c(n_size,(noccas-1),4)),
       z = zi)
}

parameters = c("phi", "p")

cjs_group_taille_trt_p <- function() {
  #likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(phi[size_group[i,t-1],t-1, Treatment[i]]*z[i,t-1])
      y[i,t] ~ dbern(p[size_group[i,t],t-1, Lake[i]]*z[i,t])
    }
  }
  #prior
  for (gs in 1:n_size){
    for (t in 1:(noccas-1)){
      for (tr in 1:4){
      phi[gs,t,tr] ~ dunif(0,1)
      }
      for (l in 1:16){
        p[gs,t,l] ~ dunif(0,1)
      }
    }
  }
}

CJS_group_taille_trt_p <- jags.parallel(data = jags.data,
                                       inits = inits,
                                       parameters.to.save = parameters,
                                       model.file = cjs_group_taille_trt_p,
                                       n.chains = 4,
                                       n.iter = ni)

save(CJS_group_taille_trt_p, file = "R/object/CJS_group_taille_trt_p.RData")


######################################## Model taille discrete


inits <- function(){
  list(p= array(rep(runif(1,0.1,0.9), n_size*(noccas-1)*16), c(n_size,(noccas-1),16)), 
       phi = array(rep(runif(1,0.1,0.9), n_size*(noccas-1)*4), c(n_size,(noccas-1),4)),
       z = zi)
}

parameters = c("phi", "p")

cjs_group_taille_trcorrect_p <- function() {
  #likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(phi[size_group[i,t-1],t-1, Treatment[i]]*z[i,t-1])
      y[i,t] ~ dbern(p[size_group[i,t],t-1, Lake[i]]*z[i,t])
    }
  }
  #prior and constraints
  for (gs in 1:n_size){
    for (t in 1:(noccas-1)){
      for (l in 1:16){
        p[gs,t,l] ~ dunif(0,1)
      }
    }
  }
  for (gs in 1:n_size){
    for (t in 1:2){
      phi[gs,t,1] ~ dunif(0,1)
      phi[gs,t,3] ~ dunif(0,1)
      phi[gs,t,2] <- phi[gs,t,1]
      phi[gs,t,4] <- phi[gs,t,3]
    }
    for (t in 3:(noccas-1)){
      for (tr in 1:4){
        phi[gs,t,tr] ~ dunif(0,1)
      }
    }
  }
}

CJS_group_taille_trcorrect_p <- jags.parallel(data = jags.data,
                                        inits = inits,
                                        parameters.to.save = parameters,
                                        model.file = cjs_group_taille_trcorrect_p,
                                        n.chains = 4,
                                        n.iter = ni)

save(CJS_group_taille_trcorrect_p, file = "R/object/CJS_group_taille_trcorrect_p.RData")

###################################### 

######################################## Model taille discrete


inits <- function(){
  list(p= runif(n_size), 
       phi = array(rep(runif(1,0.1,0.9), n_size*(noccas-1)*4), c(n_size,(noccas-1),4)),
       z = zi)
}

parameters = c("phi", "p")

cjs_group_taille_trtcorrect <- function() {
  #likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(phi[size_group[i,t-1],t-1, Treatment[i]]*z[i,t-1])
      y[i,t] ~ dbern(p[size_group[i,t]]*z[i,t])
    }
  }
  #prior and constraints
  for (gs in 1:n_size){
        p[gs] ~ dunif(0,1)
  }
  for (gs in 1:n_size){
    for (t in 1:2){
      phi[gs,t,1] ~ dunif(0,1)
      phi[gs,t,3] ~ dunif(0,1)
      phi[gs,t,2] <- phi[gs,t,1]
      phi[gs,t,4] <- phi[gs,t,3]
    }
    for (t in 3:(noccas-1)){
      for (tr in 1:4){
        phi[gs,t,tr] ~ dunif(0,1)
      }
    }
  }
}

CJS_group_taille_trtcorrect <- jags.parallel(data = jags.data,
                                              inits = inits,
                                              parameters.to.save = parameters,
                                              model.file = cjs_group_taille_trtcorrect,
                                              n.chains = 4,
                                              n.iter = ni)

save(CJS_group_taille_trtcorrect, file = "R/object/CJS_group_taille_trtcorrect.RData")

######################################## Model taille discrete


inits <- function(){
  list(p= runif(n_size), 
       phi = array(rep(runif(1,0.1,0.9), n_size*4), c(n_size,4)),
       z = zi)
}

parameters = c("phi", "p","N")

cjs_group_taille_trcorrect <- function() {
  #likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(phi_aux[size_group[i,t-1],t-1, Treatment[i]]*z[i,t-1])
      y[i,t] ~ dbern(p[size_group[i,t]]*z[i,t])
    }
  }
  #prior and constraints
  for (gs in 1:n_size){
    p[gs] ~ dunif(0,1)
  }
  for (gs in 1:n_size){
    for (tr in 1:4){
      phi[gs, tr] ~ dunif(0,1)
    }
    for (t in 1:2){
      phi_aux[gs,t,1] <- phi[gs,1]
      phi_aux[gs,t,3] <- phi[gs,1]
      
      phi_aux[gs,t,2] <- phi[gs,1]
      phi_aux[gs,t,4] <- phi[gs,3]
    }
    for (t in 3:(noccas-1)){
      for (tr in 1:4){
        phi_aux[gs,t,tr] <- phi[gs,tr]
      }
    }
  }
}

CJS_group_taille_trcorrect <- jags.parallel(data = jags.data,
                                            inits = inits,
                                            parameters.to.save = parameters,
                                            model.file = cjs_group_taille_trcorrect,
                                            n.chains = 4,
                                            n.iter = ni)

save(CJS_group_taille_trcorrect, file = "R/object/CJS_group_taille_trcorrect.RData")
