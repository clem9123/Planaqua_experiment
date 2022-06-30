#-----------------------------------------------------------------------------------------------
################################################################################################
# Model Treatment + proba_capt+ tau général sur tous les lacs
################################################################################################
#-----------------------------------------------------------------------------------------------


old = Sys.time()

model <- function(){
  
  # -------------------------------------------------
  # Parameters:
  # phi1 : survival for state 1
  # phi2 : survival for state 2
  # phi3 : survival for state 3
  # p1 : capture probability for state 1
  # p2 : capture probability for state 2
  # p3 : capture probability for state 3
  # psi12 : probability of growth from state 1 to 2
  # psi23 : probability of growth from state 2 to 3
  # error : probability of asserting a false state when captured
  # -------------------------------------------------
  # States (z):
  # 1 = alive and size below 120mm
  # 2 = alive and size between 120mm and 180mm
  # 3 = alive and size above 180mm
  # 4 = dead
  # Observations (y):  
  # 1 = detected and measured below 120mm
  # 2 = detected and measured between 120mm and 180mm
  # 3 = detected and measured above 180mm
  # 4 = not detected
  # -------------------------------------------------
  
  
  # ------
  # PRIORS
  # ------
  
  # SURVIVAL AND GROWTH
  # -------------------
  
  #survival and growth for each treatment and experimental time
  for (e in 1:2){
    for (tr in 1:4){
      phi1[tr,e] ~ dunif(0,1)
      phi2[tr,e] ~ dunif(0,1)
      phi3[tr,e] ~ dunif(0,1)
      
      psi12[tr,e] ~ dunif(0,1)
      psi23[tr,e] ~ dunif(0,1)
    }
  }
  
  for (l in 1:16) {
    lake[l]~dnorm(0,1/(sigma*sigma))
  }
  
  sigma ~ dunif(0,100)
  
  #survival and growth for each time and treatment 
  # (derived from survival and growth above, only needed for notation)
  for (t in 1:2){ #experiment time 1
    for (l in 1:16){
      logit(phi1_aux[l,t]) <- logit(phi1[LT[l],1]) + lake[l]
      logit(phi2_aux[l,t]) <- logit(phi2[LT[l],1]) + lake[l]
      logit(phi3_aux[l,t]) <- logit(phi3[LT[l],1]) + lake[l]
      
      logit(psi12_aux[l,t]) <- logit(psi12[LT[l],1]) + lake[l]
      logit(psi23_aux[l,t]) <- logit(psi23[LT[l],1]) + lake[l]
    }
  }
  
  for (t in 3:(noccas-1)){
    for (l in 1:16){ # experiment time 2
      logit(phi1_aux[l,t]) <- logit(phi1[LT[l],2]) + lake[l]
      logit(phi2_aux[l,t]) <- logit(phi2[LT[l],2]) + lake[l]
      logit(phi3_aux[l,t]) <- logit(phi3[LT[l],2]) + lake[l]
      
      logit(psi12_aux[l,t]) <- logit(psi12[LT[l],2]) + lake[l]
      logit(psi23_aux[l,t]) <- logit(psi23[LT[l],2]) + lake[l]
    }
  }
  
  # CAPTURE PROBABILITY
  # -------------------
  
  # mean capture probability
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  # variability of capture probability for each event
  for (l in 1:16){
    for (t in 1:(noccas-1)){
      epsilon[l,t] ~ dunif(-20,20)
    }
  }
  
  # corrected capture probability for each event
  for (l in 1:16){
    for(t in 1:(noccas-1)){
      logit(p1_aux[l,t]) <- logit(p1) + epsilon[l,t]
      logit(p2_aux[l,t]) <- logit(p2) + epsilon[l,t]
      logit(p3_aux[l,t]) <- logit(p3) + epsilon[l,t]
    }
  }
  
  # PROBABILITY OF ERROR WHEN INPUT
  # -------------------------------
  
  error ~ dunif(0,1)
  
  # TRANSITION PROBABILITY
  # ----------------------
  
  # probabilities of real state z(t+1) given real state z(t)
  for (t in 1:(noccas-1)){
    for (l in 1:16) {
      gamma[1,1,l,t] <- phi1_aux[l,t] * (1-psi12_aux[l,t])   # Pr(alive and small t -> alive and small t+1)    survivalin same size class
      gamma[1,2,l,t] <- phi1_aux[l,t] * psi12_aux[l,t]       # Pr(alive and small t -> alive and medium t+1)   survival and growth
      gamma[1,3,l,t] <- 0                                      # Pr(alive and small t -> alive and large t+1)    impossible to shrink
      gamma[1,4,l,t] <- (1-phi1_aux[l,t])                     # Pr(alive and small t -> dead t+1)               mortality
      gamma[2,1,l,t] <- 0                                      # Pr(alive and medium t -> alive and small t+1)   impossible to shrink
      gamma[2,2,l,t] <- phi2_aux[l,t] * (1-psi23_aux[l,t])   # Pr(alive and medium t -> alive and medium t+1)  survival in same size class
      gamma[2,3,l,t] <- phi2_aux[l,t] * psi23_aux[l,t]       # Pr(alive and medium t -> alive and large t+1)   survival and growth
      gamma[2,4,l,t] <- (1-phi2_aux[l,t])                     # Pr(alive and medium t -> dead t+1)              mortality
      gamma[3,1,l,t] <- 0                                      # Pr(alive and large t -> alive and small t+1)    impossible to shrink
      gamma[3,2,l,t] <- 0                                      # Pr(alive and large t -> alive and medium t+1)   impossible to shrink
      gamma[3,3,l,t] <- phi3_aux[l,t]                         # Pr(alive and large t -> alive and large t+1)    survivalin same size class
      gamma[3,4,l,t] <- (1-phi3_aux[l,t])                     # Pr(alive and large t -> dead t+1)               mortality
      gamma[4,1,l,t] <- 0                                      # Pr(dead t -> alive and small t+1)               no resurection
      gamma[4,2,l,t] <- 0                                      # Pr(dead t -> alive and medium t+1)              no resurection
      gamma[4,3,l,t] <- 0                                      # Pr(dead t -> alive and large t+1)               no resurection
      gamma[4,4,l,t] <- 1                                      # Pr(dead t -> dead t+1)                          dead stay s=dead (certain)
    }
  }
  
  
  # probabilities of y(t) given z(t)
  for (l in 1:16){
    for (t in 1:(noccas-1)){
      omega[1,1,l,t] <- (p1_aux[l,t]) * (1-2 * error)           # Pr( alive and small t -> captured t)
      omega[1,2,l,t] <- (p1_aux[l,t]) * error                   # Pr( alive and small t -> captured and wrongly input medium t)
      omega[1,3,l,t] <- (p1_aux[l,t]) * error                   # Pr( alive and small t -> captured  and wrongly input t)
      omega[1,4,l,t] <- 1-(p1_aux[l,t])                         # Pr( alive and small t -> not captured t)
      omega[2,1,l,t] <- (p2_aux[l,t]) * error                   # Pr( alive and medium t -> captured and wrongly input small t) 
      omega[2,2,l,t] <- (p2_aux[l,t]) * (1-2 * error)           # Pr( alive and medium t -> captured t)
      omega[2,3,l,t] <- (p2_aux[l,t]) * error                   # Pr( alive and medium t -> captured and wrongly input large t)
      omega[2,4,l,t] <- 1-(p2_aux[l,t])                         # Pr( alive and medium t -> not captured t)
      omega[3,1,l,t] <- (p3_aux[l,t]) * error                   # Pr( alive and large t -> captured and wrongly input small t)
      omega[3,2,l,t] <- (p3_aux[l,t]) * error                   # Pr( alive and large t -> captured and wrongly input medium t)
      omega[3,3,l,t] <- (p3_aux[l,t]) * (1- 2 * error)          # Pr( alive and large t -> captured t)
      omega[3,4,l,t] <- 1-(p3_aux[l,t])                         # Pr( alive and large t -> not captured t)
      omega[4,1,l,t] <- 0                                       # Pr( dead t -> captured t)
      omega[4,2,l,t] <- 0                                       # Pr( dead t -> captured t)
      omega[4,3,l,t] <- 0                                       # Pr( dead t -> captured t)
      omega[4,4,l,t] <- 1                                       # Pr( dead t -> not captured t)
    }
  }
  
  # likelihood 
  for (i in 1:nind){ # for each individual
    ## before first capture
    for (t in 1:(f[i]-1)){
      z[i,t] <- 0 # real state unknown
      
      # tables for abundance in each size class (before capture they are concidered absent)
      N1[i,t] <- 0
      N2[i,t] <- 0
      N3[i,t] <- 0
    }
    ## at first capture
    z[i,f[i]] <- fs[i] # get first capture size
    
    # input in tables for abundance in its size class
    N1[i,f[i]] <- ifelse(z[i,f[i]] == 1, 1,0)
    N2[i,f[i]] <- ifelse(z[i,f[i]] == 2, 1,0)
    N3[i,f[i]] <- ifelse(z[i,f[i]] == 3, 1,0)
    ## after first capture
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4,Lake[i],t-1])
      # input in tables for abundance in its size class
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4,Lake[i],t-1])
    }
  }
  
  # count of abundance of each size class in each treatment and year
  for (t in 1:noccas){
    #treatment 1
    n1[t,1] <- sum(N1[Tr1,t])
    n2[t,1] <- sum(N2[Tr1,t])
    n3[t,1] <- sum(N3[Tr1,t])
    ntot[t,1] <- n1[t,1] + n2[t,1] + n3[t,1]
    
    #treatment 2
    n1[t,2] <- sum(N1[Tr2,t])
    n2[t,2] <- sum(N2[Tr2,t])
    n3[t,2] <- sum(N3[Tr2,t])
    ntot[t,2] <- n1[t,2] + n2[t,2] + n3[t,2]
    
    #treatment 3
    n1[t,3] <- sum(N1[Tr3,t])
    n2[t,3] <- sum(N2[Tr3,t])
    n3[t,3] <- sum(N3[Tr3,t])
    ntot[t,3] <- n1[t,3] + n2[t,3] + n3[t,3]
    
    #treatment 4
    n1[t,4] <- sum(N1[Tr4,t])
    n2[t,4] <- sum(N2[Tr4,t])
    n3[t,4] <- sum(N3[Tr4,t])
    ntot[t,4] <- n1[t,4] + n2[t,4] + n3[t,4]
  }
}


# -------------------
# MODEL SETUP AND RUN
# -------------------

#source("R/setup.R")

inits = function(){
  list(phi1 = matrix(ncol = 2, runif(8,0,1)),phi2 = matrix(ncol = 2, runif(8,0,1)),phi3 = matrix(ncol = 2, runif(8,0,1)),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = matrix(ncol = 2, runif(8,0,0.5)),psi23 = matrix(ncol = 2, runif(8,0,0.5)), psi13 = matrix(ncol = 2, runif(8,0,0.5)),
       error = runif(1,0,0.1), epsilon = matrix(ncol = 5, runif(80,-2,2)),
       z = zi,
       sigma = runif(1,0,100))}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "error",
               "epsilon",
               "n1","n2","n3","ntot",
               "sigma")

Model_Treatment_capture_random <- jags.parallel(data = jags.data,
                                         inits = inits,
                                         parameters.to.save = parameters,
                                         model.file = model,
                                         n.chains = 2,
                                         n.iter = 10000)

save(Model_Treatment_capture_random, file = "R/Model_Treatment_capture_random.RData" )

runtime = Sys.time() - old
print(runtime)

#-----------------------------------------------------------------------------------------------
################################################################################################
# Model Treatment + tau sur tous les lacs
################################################################################################
#-----------------------------------------------------------------------------------------------


old = Sys.time()

model <- function(){
  
  # -------------------------------------------------
  # Parameters:
  # phi1 : survival for state 1
  # phi2 : survival for state 2
  # phi3 : survival for state 3
  # p1 : capture probability for state 1
  # p2 : capture probability for state 2
  # p3 : capture probability for state 3
  # psi12 : probability of growth from state 1 to 2
  # psi23 : probability of growth from state 2 to 3
  # error : probability of asserting a false state when captured
  # -------------------------------------------------
  # States (z):
  # 1 = alive and size below 120mm
  # 2 = alive and size between 120mm and 180mm
  # 3 = alive and size above 180mm
  # 4 = dead
  # Observations (y):  
  # 1 = detected and measured below 120mm
  # 2 = detected and measured between 120mm and 180mm
  # 3 = detected and measured above 180mm
  # 4 = not detected
  # -------------------------------------------------
  
  
  # ------
  # PRIORS
  # ------
  
  # SURVIVAL AND GROWTH
  # -------------------
  
  #survival and growth for each treatment and experimental time
  for (e in 1:2){
    for (tr in 1:4){
      phi1[tr,e] ~ dunif(0,1)
      phi2[tr,e] ~ dunif(0,1)
      phi3[tr,e] ~ dunif(0,1)
      
      psi12[tr,e] ~ dunif(0,1)
      psi23[tr,e] ~ dunif(0,1)
    }
  }
  
  for (l in 1:16) {
    lake[l]~dnorm(0,1/(sigma*sigma))
  }
  
  sigma ~ dunif(0,100)
  
  #survival and growth for each time and treatment 
  # (derived from survival and growth above, only needed for notation)
  for (t in 1:2){ #experiment time 1
    for (l in 1:16){
      logit(phi1_aux[l,t]) <- logit(phi1[LT[l],1]) + lake[l]
      logit(phi2_aux[l,t]) <- logit(phi2[LT[l],1]) + lake[l]
      logit(phi3_aux[l,t]) <- logit(phi3[LT[l],1]) + lake[l]
      
      logit(psi12_aux[l,t]) <- logit(psi12[LT[l],1]) + lake[l]
      logit(psi23_aux[l,t]) <- logit(psi23[LT[l],1]) + lake[l]
    }
  }
  
  for (t in 3:(noccas-1)){
    for (l in 1:16){ # experiment time 2
      logit(phi1_aux[l,t]) <- logit(phi1[LT[l],2]) + lake[l]
      logit(phi2_aux[l,t]) <- logit(phi2[LT[l],2]) + lake[l]
      logit(phi3_aux[l,t]) <- logit(phi3[LT[l],2]) + lake[l]
      
      logit(psi12_aux[l,t]) <- logit(psi12[LT[l],2]) + lake[l]
      logit(psi23_aux[l,t]) <- logit(psi23[LT[l],2]) + lake[l]
    }
  }
  
  # CAPTURE PROBABILITY
  # -------------------
  
  # mean capture probability
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  # corrected capture probability for each event
  for (l in 1:16){
    for(t in 1:(noccas-1)){
      logit(p1_aux[l,t]) <- logit(p1)
      logit(p2_aux[l,t]) <- logit(p2)
      logit(p3_aux[l,t]) <- logit(p3)
    }
  }
  
  # PROBABILITY OF ERROR WHEN INPUT
  # -------------------------------
  
  error ~ dunif(0,1)
  
  # TRANSITION PROBABILITY
  # ----------------------
  
  # probabilities of real state z(t+1) given real state z(t)
  for (t in 1:(noccas-1)){
    for (l in 1:16) {
      gamma[1,1,l,t] <- phi1_aux[l,t] * (1-psi12_aux[l,t])   # Pr(alive and small t -> alive and small t+1)    survivalin same size class
      gamma[1,2,l,t] <- phi1_aux[l,t] * psi12_aux[l,t]       # Pr(alive and small t -> alive and medium t+1)   survival and growth
      gamma[1,3,l,t] <- 0                                      # Pr(alive and small t -> alive and large t+1)    impossible to shrink
      gamma[1,4,l,t] <- (1-phi1_aux[l,t])                     # Pr(alive and small t -> dead t+1)               mortality
      gamma[2,1,l,t] <- 0                                      # Pr(alive and medium t -> alive and small t+1)   impossible to shrink
      gamma[2,2,l,t] <- phi2_aux[l,t] * (1-psi23_aux[l,t])   # Pr(alive and medium t -> alive and medium t+1)  survival in same size class
      gamma[2,3,l,t] <- phi2_aux[l,t] * psi23_aux[l,t]       # Pr(alive and medium t -> alive and large t+1)   survival and growth
      gamma[2,4,l,t] <- (1-phi2_aux[l,t])                     # Pr(alive and medium t -> dead t+1)              mortality
      gamma[3,1,l,t] <- 0                                      # Pr(alive and large t -> alive and small t+1)    impossible to shrink
      gamma[3,2,l,t] <- 0                                      # Pr(alive and large t -> alive and medium t+1)   impossible to shrink
      gamma[3,3,l,t] <- phi3_aux[l,t]                         # Pr(alive and large t -> alive and large t+1)    survivalin same size class
      gamma[3,4,l,t] <- (1-phi3_aux[l,t])                     # Pr(alive and large t -> dead t+1)               mortality
      gamma[4,1,l,t] <- 0                                      # Pr(dead t -> alive and small t+1)               no resurection
      gamma[4,2,l,t] <- 0                                      # Pr(dead t -> alive and medium t+1)              no resurection
      gamma[4,3,l,t] <- 0                                      # Pr(dead t -> alive and large t+1)               no resurection
      gamma[4,4,l,t] <- 1                                      # Pr(dead t -> dead t+1)                          dead stay s=dead (certain)
    }
  }
  
  
  # probabilities of y(t) given z(t)
  for (l in 1:16){
    for (t in 1:(noccas-1)){
      omega[1,1,l,t] <- (p1_aux[l,t]) * (1-2 * error)           # Pr( alive and small t -> captured t)
      omega[1,2,l,t] <- (p1_aux[l,t]) * error                   # Pr( alive and small t -> captured and wrongly input medium t)
      omega[1,3,l,t] <- (p1_aux[l,t]) * error                   # Pr( alive and small t -> captured  and wrongly input t)
      omega[1,4,l,t] <- 1-(p1_aux[l,t])                         # Pr( alive and small t -> not captured t)
      omega[2,1,l,t] <- (p2_aux[l,t]) * error                   # Pr( alive and medium t -> captured and wrongly input small t) 
      omega[2,2,l,t] <- (p2_aux[l,t]) * (1-2 * error)           # Pr( alive and medium t -> captured t)
      omega[2,3,l,t] <- (p2_aux[l,t]) * error                   # Pr( alive and medium t -> captured and wrongly input large t)
      omega[2,4,l,t] <- 1-(p2_aux[l,t])                         # Pr( alive and medium t -> not captured t)
      omega[3,1,l,t] <- (p3_aux[l,t]) * error                   # Pr( alive and large t -> captured and wrongly input small t)
      omega[3,2,l,t] <- (p3_aux[l,t]) * error                   # Pr( alive and large t -> captured and wrongly input medium t)
      omega[3,3,l,t] <- (p3_aux[l,t]) * (1- 2 * error)          # Pr( alive and large t -> captured t)
      omega[3,4,l,t] <- 1-(p3_aux[l,t])                         # Pr( alive and large t -> not captured t)
      omega[4,1,l,t] <- 0                                       # Pr( dead t -> captured t)
      omega[4,2,l,t] <- 0                                       # Pr( dead t -> captured t)
      omega[4,3,l,t] <- 0                                       # Pr( dead t -> captured t)
      omega[4,4,l,t] <- 1                                       # Pr( dead t -> not captured t)
    }
  }
  
  # likelihood 
  for (i in 1:nind){ # for each individual
    ## before first capture
    for (t in 1:(f[i]-1)){
      z[i,t] <- 0 # real state unknown
      
      # tables for abundance in each size class (before capture they are concidered absent)
      N1[i,t] <- 0
      N2[i,t] <- 0
      N3[i,t] <- 0
    }
    ## at first capture
    z[i,f[i]] <- fs[i] # get first capture size
    
    # input in tables for abundance in its size class
    N1[i,f[i]] <- ifelse(z[i,f[i]] == 1, 1,0)
    N2[i,f[i]] <- ifelse(z[i,f[i]] == 2, 1,0)
    N3[i,f[i]] <- ifelse(z[i,f[i]] == 3, 1,0)
    ## after first capture
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4,Lake[i],t-1])
      # input in tables for abundance in its size class
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4,Lake[i],t-1])
    }
  }
  
  # count of abundance of each size class in each treatment and year
  for (t in 1:noccas){
    #treatment 1
    n1[t,1] <- sum(N1[Tr1,t])
    n2[t,1] <- sum(N2[Tr1,t])
    n3[t,1] <- sum(N3[Tr1,t])
    ntot[t,1] <- n1[t,1] + n2[t,1] + n3[t,1]
    
    #treatment 2
    n1[t,2] <- sum(N1[Tr2,t])
    n2[t,2] <- sum(N2[Tr2,t])
    n3[t,2] <- sum(N3[Tr2,t])
    ntot[t,2] <- n1[t,2] + n2[t,2] + n3[t,2]
    
    #treatment 3
    n1[t,3] <- sum(N1[Tr3,t])
    n2[t,3] <- sum(N2[Tr3,t])
    n3[t,3] <- sum(N3[Tr3,t])
    ntot[t,3] <- n1[t,3] + n2[t,3] + n3[t,3]
    
    #treatment 4
    n1[t,4] <- sum(N1[Tr4,t])
    n2[t,4] <- sum(N2[Tr4,t])
    n3[t,4] <- sum(N3[Tr4,t])
    ntot[t,4] <- n1[t,4] + n2[t,4] + n3[t,4]
  }
}


# -------------------
# MODEL SETUP AND RUN
# -------------------

#source("R/setup.R")

inits = function(){
  list(phi1 = matrix(ncol = 2, runif(8,0,1)),phi2 = matrix(ncol = 2, runif(8,0,1)),phi3 = matrix(ncol = 2, runif(8,0,1)),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = matrix(ncol = 2, runif(8,0,0.5)),psi23 = matrix(ncol = 2, runif(8,0,0.5)), psi13 = matrix(ncol = 2, runif(8,0,0.5)),
       error = runif(1,0,0.1), #epsilon = matrix(ncol = 5, runif(80,-2,2)),
       z = zi,
       sigma = runif(1,0,100))}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "error",
               #"epsilon",
               "n1","n2","n3","ntot",
               "sigma")

Model_Treatment_random <- jags.parallel(data = jags.data,
                                                inits = inits,
                                                parameters.to.save = parameters,
                                                model.file = model,
                                                n.chains = 2,
                                                n.iter = ni)

save(Model_Treatment_random, file = "R/Model_Treatment_random.RData" )

runtime = Sys.time() - old
print(runtime)

#-----------------------------------------------------------------------------------------------
################################################################################################
# Model Treatment + proba_capt+ tau individuel par traitement
################################################################################################
#-----------------------------------------------------------------------------------------------


old = Sys.time()

model <- function(){
  
  # -------------------------------------------------
  # Parameters:
  # phi1 : survival for state 1
  # phi2 : survival for state 2
  # phi3 : survival for state 3
  # p1 : capture probability for state 1
  # p2 : capture probability for state 2
  # p3 : capture probability for state 3
  # psi12 : probability of growth from state 1 to 2
  # psi23 : probability of growth from state 2 to 3
  # error : probability of asserting a false state when captured
  # -------------------------------------------------
  # States (z):
  # 1 = alive and size below 120mm
  # 2 = alive and size between 120mm and 180mm
  # 3 = alive and size above 180mm
  # 4 = dead
  # Observations (y):  
  # 1 = detected and measured below 120mm
  # 2 = detected and measured between 120mm and 180mm
  # 3 = detected and measured above 180mm
  # 4 = not detected
  # -------------------------------------------------
  
  
  # ------
  # PRIORS
  # ------
  
  # SURVIVAL AND GROWTH
  # -------------------
  
  #survival and growth for each treatment and experimental time
  for (e in 1:2){
    for (tr in 1:4){
      phi1[tr,e] ~ dunif(0,1)
      phi2[tr,e] ~ dunif(0,1)
      phi3[tr,e] ~ dunif(0,1)
      
      psi12[tr,e] ~ dunif(0,1)
      psi23[tr,e] ~ dunif(0,1)
    }
  }
  
  for (l in 1:16) {
    lake[l]~dnorm(0,1/(sigma[LT[l]]*sigma[LT[l]]))
  }
  
  for (tr in 1:4){
    sigma[tr] ~ dunif(0,100)
  }
  
  
  #survival and growth for each time and treatment 
  # (derived from survival and growth above, only needed for notation)
  for (t in 1:2){ #experiment time 1
    for (l in 1:16){
      logit(phi1_aux[l,t]) <- logit(phi1[LT[l],1]) + lake[l]
      logit(phi2_aux[l,t]) <- logit(phi2[LT[l],1]) + lake[l]
      logit(phi3_aux[l,t]) <- logit(phi3[LT[l],1]) + lake[l]
      
      logit(psi12_aux[l,t]) <- logit(psi12[LT[l],1]) + lake[l]
      logit(psi23_aux[l,t]) <- logit(psi23[LT[l],1]) + lake[l]
    }
  }
  
  for (t in 3:(noccas-1)){
    for (l in 1:16){ # experiment time 2
      logit(phi1_aux[l,t]) <- logit(phi1[LT[l],2]) + lake[l]
      logit(phi2_aux[l,t]) <- logit(phi2[LT[l],2]) + lake[l]
      logit(phi3_aux[l,t]) <- logit(phi3[LT[l],2]) + lake[l]
      
      logit(psi12_aux[l,t]) <- logit(psi12[LT[l],2]) + lake[l]
      logit(psi23_aux[l,t]) <- logit(psi23[LT[l],2]) + lake[l]
    }
  }
  
  # CAPTURE PROBABILITY
  # -------------------
  
  # mean capture probability
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  # variability of capture probability for each event
  for (l in 1:16){
    for (t in 1:(noccas-1)){
      epsilon[l,t] ~ dunif(-20,20)
    }
  }
  
  # corrected capture probability for each event
  for (l in 1:16){
    for(t in 1:(noccas-1)){
      logit(p1_aux[l,t]) <- logit(p1) + epsilon[l,t]
      logit(p2_aux[l,t]) <- logit(p2) + epsilon[l,t]
      logit(p3_aux[l,t]) <- logit(p3) + epsilon[l,t]
    }
  }
  
  # PROBABILITY OF ERROR WHEN INPUT
  # -------------------------------
  
  error ~ dunif(0,1)
  
  # TRANSITION PROBABILITY
  # ----------------------
  
  # probabilities of real state z(t+1) given real state z(t)
  for (t in 1:(noccas-1)){
    for (l in 1:16) {
      gamma[1,1,l,t] <- phi1_aux[l,t] * (1-psi12_aux[l,t])   # Pr(alive and small t -> alive and small t+1)    survivalin same size class
      gamma[1,2,l,t] <- phi1_aux[l,t] * psi12_aux[l,t]       # Pr(alive and small t -> alive and medium t+1)   survival and growth
      gamma[1,3,l,t] <- 0                                      # Pr(alive and small t -> alive and large t+1)    impossible to shrink
      gamma[1,4,l,t] <- (1-phi1_aux[l,t])                     # Pr(alive and small t -> dead t+1)               mortality
      gamma[2,1,l,t] <- 0                                      # Pr(alive and medium t -> alive and small t+1)   impossible to shrink
      gamma[2,2,l,t] <- phi2_aux[l,t] * (1-psi23_aux[l,t])   # Pr(alive and medium t -> alive and medium t+1)  survival in same size class
      gamma[2,3,l,t] <- phi2_aux[l,t] * psi23_aux[l,t]       # Pr(alive and medium t -> alive and large t+1)   survival and growth
      gamma[2,4,l,t] <- (1-phi2_aux[l,t])                     # Pr(alive and medium t -> dead t+1)              mortality
      gamma[3,1,l,t] <- 0                                      # Pr(alive and large t -> alive and small t+1)    impossible to shrink
      gamma[3,2,l,t] <- 0                                      # Pr(alive and large t -> alive and medium t+1)   impossible to shrink
      gamma[3,3,l,t] <- phi3_aux[l,t]                         # Pr(alive and large t -> alive and large t+1)    survivalin same size class
      gamma[3,4,l,t] <- (1-phi3_aux[l,t])                     # Pr(alive and large t -> dead t+1)               mortality
      gamma[4,1,l,t] <- 0                                      # Pr(dead t -> alive and small t+1)               no resurection
      gamma[4,2,l,t] <- 0                                      # Pr(dead t -> alive and medium t+1)              no resurection
      gamma[4,3,l,t] <- 0                                      # Pr(dead t -> alive and large t+1)               no resurection
      gamma[4,4,l,t] <- 1                                      # Pr(dead t -> dead t+1)                          dead stay s=dead (certain)
    }
  }
  
  
  # probabilities of y(t) given z(t)
  for (l in 1:16){
    for (t in 1:(noccas-1)){
      omega[1,1,l,t] <- (p1_aux[l,t]) * (1-2 * error)           # Pr( alive and small t -> captured t)
      omega[1,2,l,t] <- (p1_aux[l,t]) * error                   # Pr( alive and small t -> captured and wrongly input medium t)
      omega[1,3,l,t] <- (p1_aux[l,t]) * error                   # Pr( alive and small t -> captured  and wrongly input t)
      omega[1,4,l,t] <- 1-(p1_aux[l,t])                         # Pr( alive and small t -> not captured t)
      omega[2,1,l,t] <- (p2_aux[l,t]) * error                   # Pr( alive and medium t -> captured and wrongly input small t) 
      omega[2,2,l,t] <- (p2_aux[l,t]) * (1-2 * error)           # Pr( alive and medium t -> captured t)
      omega[2,3,l,t] <- (p2_aux[l,t]) * error                   # Pr( alive and medium t -> captured and wrongly input large t)
      omega[2,4,l,t] <- 1-(p2_aux[l,t])                         # Pr( alive and medium t -> not captured t)
      omega[3,1,l,t] <- (p3_aux[l,t]) * error                   # Pr( alive and large t -> captured and wrongly input small t)
      omega[3,2,l,t] <- (p3_aux[l,t]) * error                   # Pr( alive and large t -> captured and wrongly input medium t)
      omega[3,3,l,t] <- (p3_aux[l,t]) * (1- 2 * error)          # Pr( alive and large t -> captured t)
      omega[3,4,l,t] <- 1-(p3_aux[l,t])                         # Pr( alive and large t -> not captured t)
      omega[4,1,l,t] <- 0                                       # Pr( dead t -> captured t)
      omega[4,2,l,t] <- 0                                       # Pr( dead t -> captured t)
      omega[4,3,l,t] <- 0                                       # Pr( dead t -> captured t)
      omega[4,4,l,t] <- 1                                       # Pr( dead t -> not captured t)
    }
  }
  
  # likelihood 
  for (i in 1:nind){ # for each individual
    ## before first capture
    for (t in 1:(f[i]-1)){
      z[i,t] <- 0 # real state unknown
      
      # tables for abundance in each size class (before capture they are concidered absent)
      N1[i,t] <- 0
      N2[i,t] <- 0
      N3[i,t] <- 0
    }
    ## at first capture
    z[i,f[i]] <- fs[i] # get first capture size
    
    # input in tables for abundance in its size class
    N1[i,f[i]] <- ifelse(z[i,f[i]] == 1, 1,0)
    N2[i,f[i]] <- ifelse(z[i,f[i]] == 2, 1,0)
    N3[i,f[i]] <- ifelse(z[i,f[i]] == 3, 1,0)
    ## after first capture
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4,Lake[i],t-1])
      # input in tables for abundance in its size class
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4,Lake[i],t-1])
    }
  }
  
  # count of abundance of each size class in each treatment and year
  for (t in 1:noccas){
    #treatment 1
    n1[t,1] <- sum(N1[Tr1,t])
    n2[t,1] <- sum(N2[Tr1,t])
    n3[t,1] <- sum(N3[Tr1,t])
    ntot[t,1] <- n1[t,1] + n2[t,1] + n3[t,1]
    
    #treatment 2
    n1[t,2] <- sum(N1[Tr2,t])
    n2[t,2] <- sum(N2[Tr2,t])
    n3[t,2] <- sum(N3[Tr2,t])
    ntot[t,2] <- n1[t,2] + n2[t,2] + n3[t,2]
    
    #treatment 3
    n1[t,3] <- sum(N1[Tr3,t])
    n2[t,3] <- sum(N2[Tr3,t])
    n3[t,3] <- sum(N3[Tr3,t])
    ntot[t,3] <- n1[t,3] + n2[t,3] + n3[t,3]
    
    #treatment 4
    n1[t,4] <- sum(N1[Tr4,t])
    n2[t,4] <- sum(N2[Tr4,t])
    n3[t,4] <- sum(N3[Tr4,t])
    ntot[t,4] <- n1[t,4] + n2[t,4] + n3[t,4]
  }
}


# -------------------
# MODEL SETUP AND RUN
# -------------------

#source("R/setup.R")

inits = function(){
  list(phi1 = matrix(ncol = 2, runif(8,0,1)),phi2 = matrix(ncol = 2, runif(8,0,1)),phi3 = matrix(ncol = 2, runif(8,0,1)),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = matrix(ncol = 2, runif(8,0,0.5)),psi23 = matrix(ncol = 2, runif(8,0,0.5)), psi13 = matrix(ncol = 2, runif(8,0,0.5)),
       error = runif(1,0,0.1), epsilon = matrix(ncol = 5, runif(80,-2,2)),
       z = zi,
       sigma = runif(4,0,100))}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "error",
               "epsilon",
               "n1","n2","n3","ntot",
               "sigma")

Model_Treatment_capture_random <- jags.parallel(data = jags.data,
                                                inits = inits,
                                                parameters.to.save = parameters,
                                                model.file = model,
                                                n.chains = 2,
                                                n.iter = 5000)

save(Model_Treatment_capture_randomtr, file = "R/Model_Treatment_capture_randomtr.RData" )

runtime = Sys.time() - old
print(runtime)