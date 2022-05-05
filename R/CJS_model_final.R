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

# subset of data for test
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
    if (j > f[i] & s_group[i,j]==4) {zinit[i,j] <- zinit[i,j-1]}
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

# List of index of individuals living in Lake n
#I_lake = list()
#for (i in 1:16){
#  I_lake[length(I_lake)+1] <- list(which(s$Lake == i))
#}

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
                  Tr4 = which(s$Treatment == 4))

#-----------------------------------------------------------------------------------------------
################################################################################################
# Model simple
################################################################################################
#-----------------------------------------------------------------------------------------------

multievent <- function(){
  
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
  phi1 ~ dunif(0,1)
  phi2 ~ dunif(0,1)
  phi3 ~ dunif(0,1)
  
  psi12 ~ dunif(0,1)
  psi23 ~ dunif(0,1)
  psi13 ~ dunif(0,1)
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  error ~ dunif(0,1)
  
  # probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phi1 * (1-psi12-psi13)
  gamma[1,2] <- phi1 * psi12
  gamma[1,3] <- phi1 * psi13
  gamma[1,4] <- (1-phi1)
  gamma[2,1] <- 0
  gamma[2,2] <- phi2 * (1-psi23)
  gamma[2,3] <- phi2 * psi23
  gamma[2,4] <- (1-phi2)
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- phi3
  gamma[3,4] <- (1-phi3)
  gamma[4,1] <- 0
  gamma[4,2] <- 0
  gamma[4,3] <- 0
  gamma[4,4] <- 1
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- p1 * (1-2 * error)
  omega[1,2] <- p1 * error
  omega[1,3] <- p1 * error
  omega[1,4] <- 1-p1
  omega[2,1] <- p2 * error
  omega[2,2] <- p2 * (1-error)
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
    z[i,f[i]] <- fs[i]
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
}

inits = function(){
  list(phi1 = runif(1,0,1),phi2 = runif(1,0,1),phi3 = runif(1,0,1),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = runif(1,0,0.5),psi23 = runif(1,0,0.5),
       error = runif(1,0,0.1),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12",
               "psi13","psi23",
               "error")

Model_multi <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = multievent,
                             n.chains = 2,
                             n.iter = ni)

save(Model_multi, file = "R/object/Model_multi.RData")

#-----------------------------------------------------------------------------------------------
################################################################################################
# Model simple with abundance
################################################################################################
#-----------------------------------------------------------------------------------------------

multievent_abundance <- function(){
  
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
  phi1 ~ dunif(0,1)
  phi2 ~ dunif(0,1)
  phi3 ~ dunif(0,1)
  
  psi12 ~ dunif(0,1)
  psi23 ~ dunif(0,1)
  psi13 ~ dunif(0,1)
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  error ~ dunif(0,1)
  
  # probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phi1 * (1-psi12-psi13)
  gamma[1,2] <- phi1 * psi12
  gamma[1,3] <- phi1 * psi13
  gamma[1,4] <- (1-phi1)
  gamma[2,1] <- 0
  gamma[2,2] <- phi2 * (1-psi23)
  gamma[2,3] <- phi2 * psi23
  gamma[2,4] <- (1-phi2)
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- phi3
  gamma[3,4] <- (1-phi3)
  gamma[4,1] <- 0
  gamma[4,2] <- 0
  gamma[4,3] <- 0
  gamma[4,4] <- 1
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- p1 * (1-2 * error)
  omega[1,2] <- p1 * error
  omega[1,3] <- p1 * error
  omega[1,4] <- 1-p1
  omega[2,1] <- p2 * error
  omega[2,2] <- p2 * (1-error)
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
    z[i,f[i]] <- fs[i]
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
    }
  }
  for (t in 1:noccas){
    n1[t] <- sum(N1[,t])
    n2[t] <- sum(N1[,t])
    n3[t] <- sum(N1[,t])
    ntot[t] <- n1[t] + n2[t] + n3[t]
  }
}

inits = function(){
  list(phi1 = runif(1,0,1),phi2 = runif(1,0,1),phi3 = runif(1,0,1),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = runif(1,0,0.5),psi23 = runif(1,0,0.5),
       error = runif(1,0,0.1),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12",
               "psi13","psi23",
               "error")

Model_multi_abundance <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = multievent_abundance,
                             n.chains = 2,
                             n.iter = ni)

save(Model_multi_abundance, file = "R/object/Model_multi_abundance.RData")

#-----------------------------------------------------------------------------------------------
################################################################################################
# Model treatment and time
################################################################################################
#-----------------------------------------------------------------------------------------------

multievent_treatment_time <- function(){
  
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
  for (t in 1:(noccas-1)){
    for (tr in 1:4){
      phi1[tr,t] ~ dunif(0,1)
      phi2[tr,t] ~ dunif(0,1)
      phi3[tr,t] ~ dunif(0,1)
      
      psi12[tr,t] ~ dunif(0,1)
      psi23[tr,t] ~ dunif(0,1)
      psi13[tr,t] ~ dunif(0,1)
    }
  }
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  error ~ dunif(0,1)
  
  # probabilities of state z(t+1) given z(t)
  for (t in 1:(noccas-1)){
    for (tr in 1:4) {
      gamma[1,1,tr,t] <- phi1[tr,t] * (1-psi12[tr,t]-psi13[tr,t])
      gamma[1,2,tr,t] <- phi1[tr,t] * psi12[tr,t]
      gamma[1,3,tr,t] <- phi1[tr,t] * psi13[tr,t]
      gamma[1,4,tr,t] <- (1-phi1[tr,t])
      gamma[2,1,tr,t] <- 0
      gamma[2,2,tr,t] <- phi2[tr,t] * (1-psi23[tr,t])
      gamma[2,3,tr,t] <- phi2[tr,t] * psi23[tr,t]
      gamma[2,4,tr,t] <- (1-phi2[tr,t])
      gamma[3,1,tr,t] <- 0
      gamma[3,2,tr,t] <- 0
      gamma[3,3,tr,t] <- phi3[tr,t]
      gamma[3,4,tr,t] <- (1-phi3[tr,t])
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
  omega[2,2] <- p2 * (1-error)
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
    z[i,f[i]] <- fs[i]
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4,Treatment[i],t-1])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
}

inits = function(){
  list(phi1 = matrix(ncol = 5, runif(20,0,1)),phi2 = matrix(ncol = 5, runif(20,0,1)),phi3 = matrix(ncol = 5, runif(20,0,1)),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = matrix(ncol = 5, runif(20,0,0.9)),psi23 =matrix(ncol = 5, runif(20,0,1)), psi13 = matrix(ncol = 5, runif(20,0,0.1)),
       error = runif(1,0,0.1),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "psi13",
               "error")

Model_multi_treatment_time <- jags.parallel(data = jags.data,
                                            inits = inits,
                                            parameters.to.save = parameters,
                                            model.file = multievent_treatment_time,
                                            n.chains = 2,
                                            n.iter = ni)

save(Model_multi_treatment_time, file = "R/object/Model_multi_treatment_time.RData")

#-----------------------------------------------------------------------------------------------
################################################################################################
# Model time with abundance
################################################################################################
#-----------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------
################################################################################################
# Model treatment corrected
################################################################################################
#-----------------------------------------------------------------------------------------------

multievent_treatment_corrected <- function(){
  
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
      psi13[tr,t] ~ dunif(0,1)
    }
  }
  
  for (t in 1:2){
    for (tr in 1:4){
      phi1_aux[tr,t] <- phi1[tr,1]
      phi2_aux[tr,t] <- phi2[tr,1]
      phi3_aux[tr,t] <- phi3[tr,1]
      
      psi12_aux[tr,t] <- psi12[tr,1]
      psi23_aux[tr,t] <- psi23[tr,1]
      psi13_aux[tr,t] <- psi13[tr,1]
    }
  }
  
  for (t in 3:(noccas-1)){
    for (tr in 1:4){
      phi1_aux[tr,t] <- phi1[tr,2]
      phi2_aux[tr,t] <- phi2[tr,2]
      phi3_aux[tr,t] <- phi3[tr,2]
      
      psi12_aux[tr,t] <- psi12[tr,2]
      psi23_aux[tr,t] <- psi23[tr,2]
      psi13_aux[tr,t] <- psi13[tr,2]
    }
  }
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  error ~ dunif(0,1)
  
  # probabilities of state z(t+1) given z(t)
  for (t in 1:(noccas-1)){
    for (tr in 1:4) {
      gamma[1,1,tr,t] <- phi1_aux[tr,t] * (1-psi12_aux[tr,t]-psi13_aux[tr,t])
      gamma[1,2,tr,t] <- phi1_aux[tr,t] * psi12_aux[tr,t]
      gamma[1,3,tr,t] <- phi1_aux[tr,t] * psi13_aux[tr,t]
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
  omega[2,2] <- p2 * (1-error)
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
    z[i,f[i]] <- fs[i]
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4,Treatment[i],t-1])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
}

inits = function(){
  list(phi1 = matrix(ncol = 2, runif(8,0,1)),phi2 = matrix(ncol = 2, runif(8,0,1)),phi3 = matrix(ncol = 2, runif(8,0,1)),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = matrix(ncol = 2, runif(8,0,0.5)),psi23 = matrix(ncol = 2, runif(8,0,0.5)), psi13 = matrix(ncol = 2, runif(8,0,0.5)),
       error = runif(1,0,0.1),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "psi13",
               "error")

Model_multi_treatment_corrected <- jags.parallel(data = jags.data,
                                                  inits = inits,
                                                  parameters.to.save = parameters,
                                                  model.file = multievent_treatment_corrected,
                                                  n.chains = 2,
                                                  n.iter = ni)
save(Model_multi_treatment_corrected, file = "R/object/Model_multi_treatment_corrected.RData" )

#-----------------------------------------------------------------------------------------------
################################################################################################
# Model corrected with abundance
################################################################################################
#-----------------------------------------------------------------------------------------------

multievent_treatment_corrected_abundance <- function(){
  
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
      psi13[tr,t] ~ dunif(0,1)
    }
  }
  
  for (t in 1:2){
    for (tr in 1:4){
      phi1_aux[tr,t] <- phi1[tr,1]
      phi2_aux[tr,t] <- phi2[tr,1]
      phi3_aux[tr,t] <- phi3[tr,1]
      
      psi12_aux[tr,t] <- psi12[tr,1]
      psi23_aux[tr,t] <- psi23[tr,1]
      psi13_aux[tr,t] <- psi13[tr,1]
    }
  }
  
  for (t in 3:(noccas-1)){
    for (tr in 1:4){
      phi1_aux[tr,t] <- phi1[tr,2]
      phi2_aux[tr,t] <- phi2[tr,2]
      phi3_aux[tr,t] <- phi3[tr,2]
      
      psi12_aux[tr,t] <- psi12[tr,2]
      psi23_aux[tr,t] <- psi23[tr,2]
      psi13_aux[tr,t] <- psi13[tr,2]
    }
  }
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  error ~ dunif(0,1)
  
  # probabilities of state z(t+1) given z(t)
  for (t in 1:(noccas-1)){
    for (tr in 1:4) {
      gamma[1,1,tr,t] <- phi1_aux[tr,t] * (1-psi12_aux[tr,t]-psi13_aux[tr,t])
      gamma[1,2,tr,t] <- phi1_aux[tr,t] * psi12_aux[tr,t]
      gamma[1,3,tr,t] <- phi1_aux[tr,t] * psi13_aux[tr,t]
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
  omega[2,2] <- p2 * (1-error)
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
    z[i,f[i]] <- fs[i]
    for (t in 1:f[i]){
      N1[i,t] <- ifelse(z[i,t] == 1, 1,0)
      N2[i,t] <- ifelse(z[i,t] == 2, 1,0)
      N3[i,t] <- ifelse(z[i,t] == 3, 1,0)
    }
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
    n2[t,1] <- sum(N1[Tr1,t])
    n3[t,1] <- sum(N1[Tr1,t])
    ntot[t,1] <- n1[t,1] + n2[t,1] + n3[t,1]
    
    n1[t,2] <- sum(N1[Tr2,t])
    n2[t,2] <- sum(N1[Tr2,t])
    n3[t,2] <- sum(N1[Tr2,t])
    ntot[t,2] <- n1[t,2] + n2[t,2] + n3[t,2]
    
    n1[t,3] <- sum(N1[Tr3,t])
    n2[t,3] <- sum(N1[Tr3,t])
    n3[t,3] <- sum(N1[Tr3,t])
    ntot[t,3] <- n1[t,3] + n2[t,3] + n3[t,3]
    
    n1[t,4] <- sum(N1[Tr4,t])
    n2[t,4] <- sum(N1[Tr4,t])
    n3[t,4] <- sum(N1[Tr4,t])
    ntot[t,4] <- n1[t,4] + n2[t,4] + n3[t,4]
  }
}

inits = function(){
  list(phi1 = matrix(ncol = 2, runif(8,0,1)),phi2 = matrix(ncol = 2, runif(8,0,1)),phi3 = matrix(ncol = 2, runif(8,0,1)),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = matrix(ncol = 2, runif(8,0,0.5)),psi23 = matrix(ncol = 2, runif(8,0,0.5)), psi13 = matrix(ncol = 2, runif(8,0,0.5)),
       error = runif(1,0,0.1),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "psi13",
               "error",
               "n1","n2","n3","ntot")

Model_multi_treatment_corrected_abundance <- jags.parallel(data = jags.data,
                                                  inits = inits,
                                                  parameters.to.save = parameters,
                                                  model.file = multievent_treatment_corrected_abundance,
                                                  n.chains = 2,
                                                  n.iter = ni)
save(Model_multi_treatment_corrected_abundance, file = "R/object/Model_multi_treatment_corrected_abundance.RData" )


#-----------------------------------------------------------------------------------------------
################################################################################################
# Model with time, lake and abundance : pas fonctionnel
################################################################################################
#-----------------------------------------------------------------------------------------------

multievent_lake_time_abundance <- function(){
  
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
  for (t in 1:(noccas-1)){
    for (tr in 1:16){
      phi1[tr,t] ~ dunif(0,1)
      phi2[tr,t] ~ dunif(0,1)
      phi3[tr,t] ~ dunif(0,1)
      
      psi12[tr,t] ~ dunif(0,1)
      psi23[tr,t] ~ dunif(0,1)
      psi13[tr,t] ~ dunif(0,1)
    }
  }
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  error ~ dunif(0,1)
  
  # probabilities of state z(t+1) given z(t)
  for (t in 1:(noccas-1)){
    for (tr in 1:16) {
      gamma[1,1,tr,t] <- phi1[tr,t] * (1-psi12[tr,t]-psi13[tr,t])
      gamma[1,2,tr,t] <- phi1[tr,t] * psi12[tr,t]
      gamma[1,3,tr,t] <- phi1[tr,t] * psi13[tr,t]
      gamma[1,4,tr,t] <- (1-phi1[tr,t])
      gamma[2,1,tr,t] <- 0
      gamma[2,2,tr,t] <- phi2[tr,t] * (1-psi23[tr,t])
      gamma[2,3,tr,t] <- phi2[tr,t] * psi23[tr,t]
      gamma[2,4,tr,t] <- (1-phi2[tr,t])
      gamma[3,1,tr,t] <- 0
      gamma[3,2,tr,t] <- 0
      gamma[3,3,tr,t] <- phi3[tr,t]
      gamma[3,4,tr,t] <- (1-phi3[tr,t])
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
  omega[2,2] <- p2 * (1-error)
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
    z[i,f[i]] <- fs[i]
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
    for (tr in 1:16){
      n1[t,tr] <- sum(N1[I_lake[[tr]],t])
      n2[t,tr] <- sum(N1[I_lake[[tr]],t])
      n3[t,tr] <- sum(N1[I_lake[[tr]],t])
      ntot[t,tr] <- n1[t,tr] + n2[t,tr] + n3[t,tr]
    }
  }
}

inits = function(){
  list(phi1 = matrix(ncol = 5, runif(80,0,1)),phi2 = matrix(ncol = 5, runif(80,0,1)),phi3 = matrix(ncol = 5, runif(80,0,1)),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = matrix(ncol = 5, runif(80,0,0.9)),psi23 =matrix(ncol = 5, runif(80,0,1)), psi13 = matrix(ncol = 5, runif(80,0,0.1)),
       error = runif(1,0,0.1),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "psi13",
               "error")

Model_multi_lake_time_abundance <- jags.parallel(data = jags.data,
                                            inits = inits,
                                            parameters.to.save = parameters,
                                            model.file = multievent_lake_time_abundance,
                                            n.chains = 2,
                                            n.iter = ni)

save(Model_multi_lake_time_abundance, file = "R/object/Model_lake_time_abundance.RData")