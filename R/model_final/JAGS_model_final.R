#-----------------------------------------------------------------------------------------------
################################################################################################
# Data
################################################################################################
#-----------------------------------------------------------------------------------------------

s <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from = Size)

tr <- BDD_a %>% filter(Year != "2022") %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
mooved <- tr[duplicated(tr$Tag_id),]$Tag_id
tr <- tr %>% filter ( !duplicated(tr$Tag_id))
tr <- tr %>% mutate (Treatment = ifelse(Tag_id %in% mooved, NA, Treatment))

la <- BDD_a %>% group_by(Tag_id, Lake) %>% summarize() %>% ungroup()
mooved <- la[duplicated(la$Tag_id),]$Tag_id
la <- la %>% filter ( !duplicated(la$Tag_id))
la <- la %>% mutate (Lake = ifelse(Tag_id %in% mooved, NA, Lake))

s <- merge(s,tr)
s <- merge(s,la)

s <- s %>% filter(!is.na(Lake))

#subset of data for test
#set.seed(30)
#s <- s[sample(1:nrow(s), 100),]

# creation s_group data frame
size_break <- c(0,120,180,300)
s_group <- data.frame(apply(s[2:7], 2, 
                            function(x) discretize(x,method = "fixed",
                                                   breaks = size_break, 
                                                   labels= c(1:(length(size_break)-1)))))
s_group [] <- apply(s_group, 2, as.numeric)

fs <- apply(s_group, 1, function(x) first(na.omit(x)))

s_group <- data.frame(apply(s_group,2, function(x) ifelse(is.na(x),4,x)))

get.first <- function(x) min(which(x!=4))
f <- apply(s_group, 1, get.first)

# creation zinit
zinit <- s_group
for (i in 1:nrow(s_group)) {
  for (j in 1:ncol(s_group)) {
    if (j > f[i] & s_group[i,j]==4) {zinit[i,j] <- ifelse(zinit[i,j-1]==1,2,zinit[i,j-1])}
  }
}
for (i in 1:nrow(s_group)) {
  for (j in 1:ncol(s_group)) {    
    if (j <= f[i]) {zinit[i,j] <- NA}
  }
}
for (i in 1:nrow(s_group)) {
  if (f[i]<5){
    for (j in (f[i]+2):ncol(zinit)) {
      if (zinit[i,j]<zinit[i,j-1]) {zinit[i,j] <- zinit[i,j-1]}
    }
  }}
zinit <- as.matrix(zinit)

# data
jags.data <- list(y = s_group,
                  zi = zinit,
                  f = f,
                  fs = fs,
                  nind = nrow(s_group),
                  noccas = ncol(s_group),
                  ni = 10000,
                  Lake = s$Lake,
                  Treatment = s$Treatment,
                  Tr1 = which(s$Treatment == 1),
                  Tr2 = which(s$Treatment == 2), 
                  Tr3 = which(s$Treatment == 3), 
                  Tr4 = which(s$Treatment == 4),
                  LT = Lake_treatment$Treatment)

#-----------------------------------------------------------------------------------------------
################################################################################################
# Model Treatment + lake AND time
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
  # 1 = alive and size between 
  # 2 = alive and size between
  # 3 = alive and size between
  # 4 = dead
  # Observations (y):  
  # 1 = detected and measured between
  # 2 = detected and measured between
  # 3 = detected and measured between
  # 4 = not detected
  # -------------------------------------------------
  
  # priors
  for (t in 1:2){
    for (tr in 1:4){
      phi1[tr,t] ~ dunif(0,1)
      phi2[tr,t] ~ dunif(0,1)
      phi3[tr,t] ~ dunif(0,1)
      
      psi12[tr,t] ~ dunif(0,1)
      psi23[tr,t] ~ dunif(0,1)
      #psi13[tr,t] ~ dunif(0,1)
    }
  }
  for (t in 1:(noccas-1)){
      mu_phi[t] ~ dnorm(0,1/(tau_phi^2))
      mu_psi[t] ~ dnorm(0,1/(tau_psi^2))
  }
  
  for (l in 1:16){
    epsilon_phi[l] ~ dnorm(0,1/(sigma_phi^2))
    epsilon_psi[l] ~ dnorm(0,1/(sigma_psi^2))
  }
  
  for (t in 1:2){
    for (l in 1:16){
      logit(phi1_aux[l,t]) <- logit(phi1[LT[l],1]) + epsilon_phi[l] + mu_phi[t]
      logit(phi2_aux[l,t]) <- logit(phi2[LT[l],1]) + epsilon_phi[l] + mu_phi[t]
      logit(phi3_aux[l,t]) <- logit(phi3[LT[l],1]) + epsilon_phi[l] + mu_phi[t]
      
      logit(psi12_aux[l,t]) <- logit(psi12[LT[l],1]) + epsilon_psi[l] + mu_psi[t]
      logit(psi23_aux[l,t]) <- logit(psi23[LT[l],1]) + epsilon_psi[l] + mu_psi[t]
    }
  }
  
  for (t in 3:(noccas-1)){
    for (l in 1:16){
      logit(phi1_aux[l,t]) <- logit(phi1[LT[l],2]) + epsilon_phi[l] + mu_phi[t]
      logit(phi2_aux[l,t]) <- logit(phi2[LT[l],2]) + epsilon_phi[l] + mu_phi[t]
      logit(phi3_aux[l,t]) <- logit(phi3[LT[l],2]) + epsilon_phi[l] + mu_phi[t]
      
      logit(psi12_aux[l,t]) <- logit(psi12[LT[l],2]) + epsilon_psi[l] + mu_psi[t]
      logit(psi23_aux[l,t]) <- logit(psi23[LT[l],2]) + epsilon_psi[l] + mu_psi[t]
    }
  }
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  error ~ dunif(0,1)
  
  sigma_phi ~ dunif(0,50)
  sigma_psi ~ dunif(0,50)
  tau_phi ~ dunif(0,50)
  tau_psi ~ dunif(0,50)
  
  # probabilities of state z(t+1) given z(t)
  for (t in 1:(noccas-1)){
    for (l in 1:16) {
      gamma[1,1,l,t] <- phi1_aux[l,t] * (1-psi12_aux[l,t])#-psi13_aux[l,t])
      gamma[1,2,l,t] <- phi1_aux[l,t] * psi12_aux[l,t]
      gamma[1,3,l,t] <- 0 #phi1_aux[l,t] * psi13_aux[l,t]
      gamma[1,4,l,t] <- (1-phi1_aux[l,t])
      gamma[2,1,l,t] <- 0
      gamma[2,2,l,t] <- phi2_aux[l,t] * (1-psi23_aux[l,t])
      gamma[2,3,l,t] <- phi2_aux[l,t] * psi23_aux[l,t]
      gamma[2,4,l,t] <- (1-phi2_aux[l,t])
      gamma[3,1,l,t] <- 0
      gamma[3,2,l,t] <- 0
      gamma[3,3,l,t] <- phi3_aux[l,t]
      gamma[3,4,l,t] <- (1-phi3_aux[l,t])
      gamma[4,1,l,t] <- 0
      gamma[4,2,l,t] <- 0
      gamma[4,3,l,t] <- 0
      gamma[4,4,l,t] <- 1
    }
  }
  
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- p1 * (1-2 * error)
  omega[1,2] <- p1 * error
  omega[1,3] <- p1 * error
  omega[1,4] <- 1-p1
  omega[2,1] <- p2 * error
  omega[2,2] <- p2 * (1-2*error)
  omega[2,3] <- p2 * error
  omega[2,4] <- 1-p2
  omega[3,1] <- p3 * error
  omega[3,2] <- p3 * error
  omega[3,3] <- p3 * (1- 2*error)
  omega[3,4] <- 1-p3
  omega[4,1] <- 0
  omega[4,2] <- 0
  omega[4,3] <- 0
  omega[4,4] <- 1
  
  # likelihood 
  for (i in 1:nind){
    # State at first capture
    
    for (t in 1:(f[i]-1)){
      z[i,t] <- 0
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
    }
    z[i,f[i]] <- fs[i]
    N1[i,f[i]] <- ifelse(z[i,f[i]] == 1, 1,0)
    N2[i,f[i]] <- ifelse(z[i,f[i]] == 2, 1,0)
    N3[i,f[i]] <- ifelse(z[i,f[i]] == 3, 1,0)
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4,Lake[i],t-1])
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
  for (t in 1:noccas){
    n1[t,1] <- sum(N1[Tr1,t])
    n2[t,1] <- sum(N2[Tr1,t])
    n3[t,1] <- sum(N3[Tr1,t])
    ntot[t,1] <- n1[t,1] + n2[t,1] + n3[t,1]
    
    n1[t,2] <- sum(N1[Tr2,t])
    n2[t,2] <- sum(N2[Tr2,t])
    n3[t,2] <- sum(N3[Tr2,t])
    ntot[t,2] <- n1[t,2] + n2[t,2] + n3[t,2]
    
    n1[t,3] <- sum(N1[Tr3,t])
    n2[t,3] <- sum(N2[Tr3,t])
    n3[t,3] <- sum(N3[Tr3,t])
    ntot[t,3] <- n1[t,3] + n2[t,3] + n3[t,3]
    
    n1[t,4] <- sum(N1[Tr4,t])
    n2[t,4] <- sum(N2[Tr4,t])
    n3[t,4] <- sum(N3[Tr4,t])
    ntot[t,4] <- n1[t,4] + n2[t,4] + n3[t,4]
  }
}

inits = function(){
  list(phi1 = matrix(ncol = 2, runif(8,0,1)),phi2 = matrix(ncol = 2, runif(8,0,1)),phi3 = matrix(ncol = 2, runif(8,0,1)),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = matrix(ncol = 2, runif(8,0,0.5)),psi23 = matrix(ncol = 2, runif(8,0,0.5)), #psi13 = matrix(ncol = 2, runif(8,0,0.5)),
       error = runif(1,0,0.1),
       sigma_phi = runif(1,0,1), sigma_psi = runif(1,0,1),
       tau_phi = runif(1,0,1), tau_psi = runif(1,0,1),
       epsilon_psi = runif(16, -1,1),mu_psi = runif(5,-1,1),
       epsilon_phi = runif(16, -1,1),mu_phi = runif(5,-1,1),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "epsilon_phi","mu_phi","sigma_phi","tau_phi",
               "epsilon_psi","mu_psi","sigma_psi","tau_psi",
               "error", 
               "n1","n2","n3","ntot")

Model_Treatment_time_and_lake <- jags.parallel(data = jags.data,
                                               inits = inits,
                                               parameters.to.save = parameters,
                                               model.file = model,
                                               n.chains = 2,
                                               n.iter = ni)
save(Model_Treatment_time_and_lake, file = "R/model_final/Model_Treatment_time_and_lake.RData" )

runtime = Sys.time() - old
print(runtime)

#-----------------------------------------------------------------------------------------------
################################################################################################
# Model Treatment + lake X time
################################################################################################
#-----------------------------------------------------------------------------------------------

old = Sys.time()
print(old)

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
  # 1 = alive and size between 
  # 2 = alive and size between
  # 3 = alive and size between
  # 4 = dead
  # Observations (y):  
  # 1 = detected and measured between
  # 2 = detected and measured between
  # 3 = detected and measured between
  # 4 = not detected
  # -------------------------------------------------
  
  # priors
  for (t in 1:2){
    for (tr in 1:4){
      phi1[tr,t] ~ dunif(0,1)
      phi2[tr,t] ~ dunif(0,1)
      phi3[tr,t] ~ dunif(0,1)
      
      psi12[tr,t] ~ dunif(0,1)
      psi23[tr,t] ~ dunif(0,1)
      #psi13[tr,t] ~ dunif(0,1)
    }
  }
  for (t in 1:(noccas-1)){
    for (l in 1:16){
      epsilon_phi[l,t] ~ dnorm(0,1/(sigma_phi^2))
      epsilon_psi[l,t] ~ dnorm(0,1/(sigma_psi^2))
    }
  }
  
  for (t in 1:2){
    for (l in 1:16){
      logit(phi1_aux[l,t]) <- logit(phi1[LT[l],1]) + epsilon_phi[l,t]
      logit(phi2_aux[l,t]) <- logit(phi2[LT[l],1]) + epsilon_phi[l,t]
      logit(phi3_aux[l,t]) <- logit(phi3[LT[l],1]) + epsilon_phi[l,t]
      
      logit(psi12_aux[l,t]) <- logit(psi12[LT[l],1]) + epsilon_psi[l,t]
      logit(psi23_aux[l,t]) <- logit(psi23[LT[l],1]) + epsilon_psi[l,t]
      #psi13_aux[l,t] ~ dnorm(psi13[tr,1], 1/(sigma^2))
    }
  }
  
  for (t in 3:(noccas-1)){
    for (l in 1:16){
      logit(phi1_aux[l,t]) <- logit(phi1[LT[l],2]) + epsilon_phi[l,t]
      logit(phi2_aux[l,t]) <- logit(phi2[LT[l],2]) + epsilon_phi[l,t]
      logit(phi3_aux[l,t]) <- logit(phi3[LT[l],2]) + epsilon_phi[l,t]
      
      logit(psi12_aux[l,t]) <- logit(psi12[LT[l],2]) + epsilon_psi[l,t]
      logit(psi23_aux[l,t]) <- logit(psi23[LT[l],2]) + epsilon_psi[l,t]
      #psi13_aux[l,t] ~ dnorm(psi13[tr,2], 1/(sigma^2))
    }
  }
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  error ~ dunif(0,1)
  
  sigma_phi ~ dunif(0,50)
  sigma_psi ~ dunif(0,50)
  
  # probabilities of state z(t+1) given z(t)
  for (t in 1:(noccas-1)){
    for (l in 1:16) {
      gamma[1,1,l,t] <- phi1_aux[l,t] * (1-psi12_aux[l,t])#-psi13_aux[l,t])
      gamma[1,2,l,t] <- phi1_aux[l,t] * psi12_aux[l,t]
      gamma[1,3,l,t] <- 0 #phi1_aux[l,t] * psi13_aux[l,t]
      gamma[1,4,l,t] <- (1-phi1_aux[l,t])
      gamma[2,1,l,t] <- 0
      gamma[2,2,l,t] <- phi2_aux[l,t] * (1-psi23_aux[l,t])
      gamma[2,3,l,t] <- phi2_aux[l,t] * psi23_aux[l,t]
      gamma[2,4,l,t] <- (1-phi2_aux[l,t])
      gamma[3,1,l,t] <- 0
      gamma[3,2,l,t] <- 0
      gamma[3,3,l,t] <- phi3_aux[l,t]
      gamma[3,4,l,t] <- (1-phi3_aux[l,t])
      gamma[4,1,l,t] <- 0
      gamma[4,2,l,t] <- 0
      gamma[4,3,l,t] <- 0
      gamma[4,4,l,t] <- 1
    }
  }
  
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- p1 * (1-2 * error)
  omega[1,2] <- p1 * error
  omega[1,3] <- p1 * error
  omega[1,4] <- 1-p1
  omega[2,1] <- p2 * error
  omega[2,2] <- p2 * (1- 2*error)
  omega[2,3] <- p2 * error
  omega[2,4] <- 1-p2
  omega[3,1] <- p3 * error
  omega[3,2] <- p3 * error
  omega[3,3] <- p3 * (1- 2*error)
  omega[3,4] <- 1-p3
  omega[4,1] <- 0
  omega[4,2] <- 0
  omega[4,3] <- 0
  omega[4,4] <- 1
  
  # likelihood 
  for (i in 1:nind){
    # State at first capture
    
    for (t in 1:(f[i]-1)){
      z[i,t] <- 0
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
    }
    z[i,f[i]] <- fs[i]
    N1[i,f[i]] <- ifelse(z[i,f[i]] == 1, 1,0)
    N2[i,f[i]] <- ifelse(z[i,f[i]] == 2, 1,0)
    N3[i,f[i]] <- ifelse(z[i,f[i]] == 3, 1,0)
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4,Lake[i],t-1])
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
  for (t in 1:noccas){
    n1[t,1] <- sum(N1[Tr1,t])
    n2[t,1] <- sum(N2[Tr1,t])
    n3[t,1] <- sum(N3[Tr1,t])
    ntot[t,1] <- n1[t,1] + n2[t,1] + n3[t,1]
    
    n1[t,2] <- sum(N1[Tr2,t])
    n2[t,2] <- sum(N2[Tr2,t])
    n3[t,2] <- sum(N3[Tr2,t])
    ntot[t,2] <- n1[t,2] + n2[t,2] + n3[t,2]
    
    n1[t,3] <- sum(N1[Tr3,t])
    n2[t,3] <- sum(N2[Tr3,t])
    n3[t,3] <- sum(N3[Tr3,t])
    ntot[t,3] <- n1[t,3] + n2[t,3] + n3[t,3]
    
    n1[t,4] <- sum(N1[Tr4,t])
    n2[t,4] <- sum(N2[Tr4,t])
    n3[t,4] <- sum(N3[Tr4,t])
    ntot[t,4] <- n1[t,4] + n2[t,4] + n3[t,4]
  }
}

inits = function(){
  list(phi1 = matrix(ncol = 2, runif(8,0,1)),phi2 = matrix(ncol = 2, runif(8,0,1)),phi3 = matrix(ncol = 2, runif(8,0,1)),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = matrix(ncol = 2, runif(8,0,0.5)),psi23 = matrix(ncol = 2, runif(8,0,0.5)), #psi13 = matrix(ncol = 2, runif(8,0,0.5)),
       error = runif(1,0,0.1), sigma_phi = runif(1,0,1),sigma_psi = runif(1,0,1),
       epsilon_phi = array(runif(80,0,1),c(16,5)), epsilon_psi = array(runif(80,0,1),c(16,5)),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "epsilon_phi","sigma_phi",
               "epsilon_psi","sigma_psi",
               "error", 
               "n1","n2","n3","ntot")

Model_Treatment_time_x_lake <- jags.parallel(data = jags.data,
                                              inits = inits,
                                              parameters.to.save = parameters,
                                              model.file = model,
                                              n.chains = 2,
                                              n.iter = ni)
save(Model_Treatment_time_x_lake, file = "R/model_final/Model_Treatment_time_x_lake.RData" )

runtime = Sys.time() - old
print(runtime)

#-----------------------------------------------------------------------------------------------
################################################################################################
# Model Treatment
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
  # 1 = alive and size between 
  # 2 = alive and size between
  # 3 = alive and size between
  # 4 = dead
  # Observations (y):  
  # 1 = detected and measured between
  # 2 = detected and measured between
  # 3 = detected and measured between
  # 4 = not detected
  # -------------------------------------------------
  
  # priors
  for (e in 1:2){
    for (tr in 1:4){
      phi1[tr,e] ~ dunif(0,1)
      phi2[tr,e] ~ dunif(0,1)
      phi3[tr,e] ~ dunif(0,1)
      
      psi12[tr,e] ~ dunif(0,1)
      psi23[tr,e] ~ dunif(0,1)
      #psi13[tr,t] ~ dunif(0,1)
    }
  }
  
  for (t in 1:2){
    for (tr in 1:4){
      phi1_aux[tr,t] <- phi1[tr,1]
      phi2_aux[tr,t] <- phi2[tr,1]
      phi3_aux[tr,t] <- phi3[tr,1]
      
      psi12_aux[tr,t] <- psi12[tr,1]
      psi23_aux[tr,t] <- psi23[tr,1]
      #psi13_aux[tr,t] <- psi13[tr,1]
    }
  }
  
  for (t in 3:(noccas-1)){
    for (tr in 1:4){
      phi1_aux[tr,t] <- phi1[tr,2]
      phi2_aux[tr,t] <- phi2[tr,2]
      phi3_aux[tr,t] <- phi3[tr,2]
      
      psi12_aux[tr,t] <- psi12[tr,2]
      psi23_aux[tr,t] <- psi23[tr,2]
      #psi13_aux[tr,t] <- psi13[tr,2]
    }
  }
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  error ~ dunif(0,1)
  
  # probabilities of state z(t+1) given z(t)
  for (t in 1:(noccas-1)){
    for (tr in 1:4) {
      gamma[1,1,tr,t] <- phi1_aux[tr,t] * (1-psi12_aux[tr,t]) #-psi13_aux[tr,t])
      gamma[1,2,tr,t] <- phi1_aux[tr,t] * psi12_aux[tr,t]
      gamma[1,3,tr,t] <- 0 #phi1_aux[tr,t] * psi13_aux[tr,t]
      gamma[1,4,tr,t] <- (1-phi1_aux[tr,t])
      gamma[2,1,tr,t] <- 0
      gamma[2,2,tr,t] <- phi2_aux[tr,t] * (1-psi23_aux[tr,t])
      gamma[2,3,tr,t] <- phi2_aux[tr,t] * psi23_aux[tr,t]
      gamma[2,4,tr,t] <- (1-phi2_aux[tr,t])
      gamma[3,1,tr,t] <- 0
      gamma[3,2,tr,t] <- 0
      gamma[3,3,tr,t] <- phi3_aux[tr,t]
      gamma[3,4,tr,t] <- (1-phi3_aux[tr,t])
      gamma[4,1,tr,t] <- 0
      gamma[4,2,tr,t] <- 0
      gamma[4,3,tr,t] <- 0
      gamma[4,4,tr,t] <- 1
    }
  }
  
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- p1 * (1-2 * error)
  omega[1,2] <- p1 * error
  omega[1,3] <- p1 * error
  omega[1,4] <- 1-p1
  omega[2,1] <- p2 * error
  omega[2,2] <- p2 * (1-2 * error)
  omega[2,3] <- p2 * error
  omega[2,4] <- 1-p2
  omega[3,1] <- p3 * error
  omega[3,2] <- p3 * error
  omega[3,3] <- p3 * (1- 2*error)
  omega[3,4] <- 1-p3
  omega[4,1] <- 0
  omega[4,2] <- 0
  omega[4,3] <- 0
  omega[4,4] <- 1
  
  # likelihood 
  for (i in 1:nind){
    # State at first capture
    
    for (t in 1:(f[i]-1)){
      z[i,t] <- 0
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
    }
    z[i,f[i]] <- fs[i]
    N1[i,f[i]] <- ifelse(z[i,f[i]] == 1, 1,0)
    N2[i,f[i]] <- ifelse(z[i,f[i]] == 2, 1,0)
    N3[i,f[i]] <- ifelse(z[i,f[i]] == 3, 1,0)
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4,Treatment[i],t-1])
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
  for (t in 1:noccas){
    n1[t,1] <- sum(N1[Tr1,t])
    n2[t,1] <- sum(N2[Tr1,t])
    n3[t,1] <- sum(N3[Tr1,t])
    ntot[t,1] <- n1[t,1] + n2[t,1] + n3[t,1]
    
    n1[t,2] <- sum(N1[Tr2,t])
    n2[t,2] <- sum(N2[Tr2,t])
    n3[t,2] <- sum(N3[Tr2,t])
    ntot[t,2] <- n1[t,2] + n2[t,2] + n3[t,2]
    
    n1[t,3] <- sum(N1[Tr3,t])
    n2[t,3] <- sum(N2[Tr3,t])
    n3[t,3] <- sum(N3[Tr3,t])
    ntot[t,3] <- n1[t,3] + n2[t,3] + n3[t,3]
    
    n1[t,4] <- sum(N1[Tr4,t])
    n2[t,4] <- sum(N2[Tr4,t])
    n3[t,4] <- sum(N3[Tr4,t])
    ntot[t,4] <- n1[t,4] + n2[t,4] + n3[t,4]
  }
}

inits = function(){
  list(phi1 = matrix(ncol = 2, runif(8,0,1)),phi2 = matrix(ncol = 2, runif(8,0,1)),phi3 = matrix(ncol = 2, runif(8,0,1)),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = matrix(ncol = 2, runif(8,0,0.5)),psi23 = matrix(ncol = 2, runif(8,0,0.5)),
       error = runif(1,0,0.1),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "error",
               "n1","n2","n3","ntot")

Model_Treatment <- jags.parallel(data = jags.data,
                                 inits = inits,
                                 parameters.to.save = parameters,
                                 model.file = model,
                                 n.chains = 2,
                                 n.iter = ni)

save(Model_Treatment, file = "R/model_final/Model_Treatment.RData" )

runtime = Sys.time() - old
print(runtime)

#-----------------------------------------------------------------------------------------------
################################################################################################
# Model Treatment + proba_capt
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
  # 1 = alive and size between 
  # 2 = alive and size between
  # 3 = alive and size between
  # 4 = dead
  # Observations (y):  
  # 1 = detected and measured between
  # 2 = detected and measured between
  # 3 = detected and measured between
  # 4 = not detected
  # -------------------------------------------------
  
  # priors
  for (e in 1:2){
    for (tr in 1:4){
      phi1[tr,e] ~ dunif(0,1)
      phi2[tr,e] ~ dunif(0,1)
      phi3[tr,e] ~ dunif(0,1)
      
      psi12[tr,e] ~ dunif(0,1)
      psi23[tr,e] ~ dunif(0,1)
      #psi13[tr,t] ~ dunif(0,1)
    }
  }
  
  for (t in 1:2){
    for (tr in 1:4){
      phi1_aux[tr,t] <- phi1[tr,1]
      phi2_aux[tr,t] <- phi2[tr,1]
      phi3_aux[tr,t] <- phi3[tr,1]
      
      psi12_aux[tr,t] <- psi12[tr,1]
      psi23_aux[tr,t] <- psi23[tr,1]
      #psi13_aux[tr,t] <- psi13[tr,1]
    }
  }
  
  for (t in 3:(noccas-1)){
    for (tr in 1:4){
      phi1_aux[tr,t] <- phi1[tr,2]
      phi2_aux[tr,t] <- phi2[tr,2]
      phi3_aux[tr,t] <- phi3[tr,2]
      
      psi12_aux[tr,t] <- psi12[tr,2]
      psi23_aux[tr,t] <- psi23[tr,2]
      #psi13_aux[tr,t] <- psi13[tr,2]
    }
  }
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  for (l in 1:16){
    for (t in 1:(noccas-1)){
      epsilon[l,t] ~ dunif(-20,20)
    }
  }
  
  for (l in 1:16){
    for(t in 1:(noccas-1)){
      logit(p1_aux[l,t]) <- logit(p1) + epsilon[l,t]
      logit(p2_aux[l,t]) <- logit(p2) + epsilon[l,t]
      logit(p3_aux[l,t]) <- logit(p3) + epsilon[l,t]
    }
  }

  error ~ dunif(0,1)
  
  # probabilities of state z(t+1) given z(t)
  for (t in 1:(noccas-1)){
    for (tr in 1:4) {
      gamma[1,1,tr,t] <- phi1_aux[tr,t] * (1-psi12_aux[tr,t]) #-psi13_aux[tr,t])
      gamma[1,2,tr,t] <- phi1_aux[tr,t] * psi12_aux[tr,t]
      gamma[1,3,tr,t] <- 0 #phi1_aux[tr,t] * psi13_aux[tr,t]
      gamma[1,4,tr,t] <- (1-phi1_aux[tr,t])
      gamma[2,1,tr,t] <- 0
      gamma[2,2,tr,t] <- phi2_aux[tr,t] * (1-psi23_aux[tr,t])
      gamma[2,3,tr,t] <- phi2_aux[tr,t] * psi23_aux[tr,t]
      gamma[2,4,tr,t] <- (1-phi2_aux[tr,t])
      gamma[3,1,tr,t] <- 0
      gamma[3,2,tr,t] <- 0
      gamma[3,3,tr,t] <- phi3_aux[tr,t]
      gamma[3,4,tr,t] <- (1-phi3_aux[tr,t])
      gamma[4,1,tr,t] <- 0
      gamma[4,2,tr,t] <- 0
      gamma[4,3,tr,t] <- 0
      gamma[4,4,tr,t] <- 1
    }
  }
  
  
  # probabilities of y(t) given z(t)
  for (l in 1:16){
    for (t in 1:(noccas-1)){
      omega[1,1,l,t] <- (p1_aux[l,t]) * (1-2 * error)
      omega[1,2,l,t] <- (p1_aux[l,t]) * error
      omega[1,3,l,t] <- (p1_aux[l,t]) * error
      omega[1,4,l,t] <- 1-(p1_aux[l,t])
      omega[2,1,l,t] <- (p2_aux[l,t]) * error
      omega[2,2,l,t] <- (p2_aux[l,t]) * (1-2 * error)
      omega[2,3,l,t] <- (p2_aux[l,t]) * error
      omega[2,4,l,t] <- 1-(p2_aux[l,t])
      omega[3,1,l,t] <- (p3_aux[l,t]) * error
      omega[3,2,l,t] <- (p3_aux[l,t]) * error
      omega[3,3,l,t] <- (p3_aux[l,t]) * (1- 2 * error)
      omega[3,4,l,t] <- 1-(p3_aux[l,t])
      omega[4,1,l,t] <- 0
      omega[4,2,l,t] <- 0
      omega[4,3,l,t] <- 0
      omega[4,4,l,t] <- 1
    }
  }
  
  # likelihood 
  for (i in 1:nind){
    # State at first capture
    
    for (t in 1:(f[i]-1)){
      z[i,t] <- 0
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
    }
    z[i,f[i]] <- fs[i]
    N1[i,f[i]] <- ifelse(z[i,f[i]] == 1, 1,0)
    N2[i,f[i]] <- ifelse(z[i,f[i]] == 2, 1,0)
    N3[i,f[i]] <- ifelse(z[i,f[i]] == 3, 1,0)
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4,Treatment[i],t-1])
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4,Lake[i],t-1])
    }
  }
  for (t in 1:noccas){
    n1[t,1] <- sum(N1[Tr1,t])
    n2[t,1] <- sum(N2[Tr1,t])
    n3[t,1] <- sum(N3[Tr1,t])
    ntot[t,1] <- n1[t,1] + n2[t,1] + n3[t,1]
    
    n1[t,2] <- sum(N1[Tr2,t])
    n2[t,2] <- sum(N2[Tr2,t])
    n3[t,2] <- sum(N3[Tr2,t])
    ntot[t,2] <- n1[t,2] + n2[t,2] + n3[t,2]
    
    n1[t,3] <- sum(N1[Tr3,t])
    n2[t,3] <- sum(N2[Tr3,t])
    n3[t,3] <- sum(N3[Tr3,t])
    ntot[t,3] <- n1[t,3] + n2[t,3] + n3[t,3]
    
    n1[t,4] <- sum(N1[Tr4,t])
    n2[t,4] <- sum(N2[Tr4,t])
    n3[t,4] <- sum(N3[Tr4,t])
    ntot[t,4] <- n1[t,4] + n2[t,4] + n3[t,4]
  }
}

inits = function(){
  list(phi1 = matrix(ncol = 2, runif(8,0,1)),phi2 = matrix(ncol = 2, runif(8,0,1)),phi3 = matrix(ncol = 2, runif(8,0,1)),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = matrix(ncol = 2, runif(8,0,0.5)),psi23 = matrix(ncol = 2, runif(8,0,0.5)), psi13 = matrix(ncol = 2, runif(8,0,0.5)),
       error = runif(1,0,0.1), epsilon = matrix(ncol = 5, runif(80,-2,2)),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "error",
               "epsilon",
               "n1","n2","n3","ntot")

Model_treatment_capture <- jags.parallel(data = jags.data,
                                 inits = inits,
                                 parameters.to.save = parameters,
                                 model.file = model,
                                 n.chains = 2,
                                 n.iter = ni)

save(Model_Treatment_capture, file = "R/model_final/Model_Treatment_capture.RData" )

runtime = Sys.time() - old
print(runtime)

