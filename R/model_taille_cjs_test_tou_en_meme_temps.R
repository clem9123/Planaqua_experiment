cjs_taille1 <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      p[i,t_bis] <- 0
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      s[i,t] <- size[i,t]
      logit(p[i,t]) <- a_p*10^-2*s[i,t] + b_p #a2_p*10^-4*(s[i,t]^2) + 
      logit(phi[i,t]) <- a_phi*10^-2*s[i,t] + b_phi # + a2_phi*10^-4*size[i,t]^2
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  a_p ~ dunif(-1000,1000)
  a_phi ~ dunif(-1000,1000)
  #a2_p ~ dunif(-1000,1000)
  #a2_phi ~ dunif(-10,0)
  b_phi ~ dunif(-1000,1000)
  b_p ~ dunif(-1000,1000)
 
  #Priors
  sd <- 12 #~ dunif(0,100)
  K <- 0.352 #~ dunif(-1000,1000)
  Linf <- 187 #~ dunif(0,500)
  #Likelihood
  for (i in 1:nind){
    for (t in (f[i]+1):noccas){
      size[i,t] ~ dnorm(size[i,t-1] + K*(Linf-size[i,t-1]), 1/(sd^2))
    } # t
  } # i
}

inits <- function(){
  list(a_p = runif(1,5,13), a_phi = runif(1,1,20),
       #a2_p = runif(1,-7,0), #a2_phi = runif(1,-7,0),
       b_p = runif(1,-10,5), b_phi = runif(1,-20,10),
       z = zi)}#,
       #sd = runif(1,0,100),
       #g =runif(1,-100,100), 
       #Linf = runif(1,0,500)
parameters = c("a2_p","a_p","a_phi","b_p","b_phi","K","sd","Linf")

CJS_taille1 <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = cjs_taille1,
                             n.chains = 3,
                             n.iter = 1000)

taille <- function(){
  #Priors
  sd ~ dunif(0,100)
  for (tr in 1:4){
    K[tr] ~ dunif(-1000,1000)
    Linf[tr] ~ dunif(0,500)
  }
  
  #Likelihood
  for (i in 1:nind){
    for (t in (f[i]+1):noccas){
      size[i,t] ~ dnorm(size[i,t-1] + K[Treatment[i]]*(Linf[Treatment[i]]-size[i,t-1]), 1/(sd^2))
    } # t
  } # i
}

inits = function(){
  list(sd = runif(1,0,100), K = runif(4,-1000,1000), Linf = runif(4, 0,500))
}

parameters = c("K","Linf","sd")

Taille <- jags.parallel(data = jags.data,
               inits = inits,
               parameters.to.save = parameters,
               model.file = taille,
               n.iter = 1000,
               n.chains = 2)


#####################################
#####################################

old <- Sys.time()

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

#set.seed(30)
#s <- s[sample(1:nrow(s), 200),]

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

# création de zinit

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


# mise en forme des données
jags.data <- list(y = s_group,
                  zi = zinit,
                  f = f,
                  fs = fs,
                  nind = nrow(s_group),
                  noccas = ncol(s_group),
                  ni = 10000,
                  Lake = s$Lake,
                  Treatment = s$Treatment) 

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

runtime1 = Sys.time() - old

multievent_treatment <- function(){
  
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
  for (tr in 1:4){
    phi1[tr] ~ dunif(0,1)
    phi2[tr] ~ dunif(0,1)
    phi3[tr] ~ dunif(0,1)
    
    psi12[tr] ~ dunif(0,1)
    psi23[tr] ~ dunif(0,1)
    psi13[tr] ~ dunif(0,1)
  }
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
    
  error ~ dunif(0,1)
  
  # probabilities of state z(t+1) given z(t)

  for (tr in 1:4) {
    gamma[1,1,tr] <- phi1[tr] * (1-psi12[tr]-psi13[tr])
    gamma[1,2,tr] <- phi1[tr] * psi12[tr]
    gamma[1,3,tr] <- phi1[tr] * psi13[tr]
    gamma[1,4,tr] <- (1-phi1[tr])
    gamma[2,1,tr] <- 0
    gamma[2,2,tr] <- phi2[tr] * (1-psi23[tr])
    gamma[2,3,tr] <- phi2[tr] * psi23[tr]
    gamma[2,4,tr] <- (1-phi2[tr])
    gamma[3,1,tr] <- 0
    gamma[3,2,tr] <- 0
    gamma[3,3,tr] <- phi3[tr]
    gamma[3,4,tr] <- (1-phi3[tr])
    gamma[4,1,tr] <- 0
    gamma[4,2,tr] <- 0
    gamma[4,3,tr] <- 0
    gamma[4,4,tr] <- 1
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
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4,Treatment[i]])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
}

inits = function(){
  list(phi1 = runif(4,0,1),phi2 = runif(4,0,1),phi3 = runif(4,0,1),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = runif(4,0,0.5),psi23 = runif(4,0,0.5), psi13 = runif(4,0,0.5),
       error = runif(1,0,0.1),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi23",
               "psi13",
               "error")

Model_multi_treatment <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = multievent_treatment,
                             n.chains = 2,
                             n.iter = ni)

runtime2 = Sys.time() - old - runtime1

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

runtime3 = Sys.time() - old - runtime2 - runtime1

multievent_treatment_time_corrected <- function(){
  
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
    for (tr in c(2,4)){
      phi1[tr,t] ~ dunif(0,1)
      phi2[tr,t] ~ dunif(0,1)
      phi3[tr,t] ~ dunif(0,1)
      
      psi12[tr,t] ~ dunif(0,1)
      psi23[tr,t] ~ dunif(0,1)
      psi13[tr,t] ~ dunif(0,1)
    }
    for (tr in c(1,3)){
      phi1[tr,t] <- phi1[tr+1,t]
      phi2[tr,t] <- phi2[tr+1,t]
      phi3[tr,t] <- phi3[tr+1,t]
      
      psi12[tr,t] <- psi12[tr+1,t]
      psi23[tr,t] <- psi23[tr+1,t]
      psi13[tr,t] <- psi13[tr+1,t]
    }
  }
  
  for (t in 3:(noccas-1)){
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

Model_multi_treatment_time_corrected <- jags.parallel(data = jags.data,
                                            inits = inits,
                                            parameters.to.save = parameters,
                                            model.file = multievent_treatment_time_corrected,
                                            n.chains = 2,
                                            n.iter = ni)
save(Model_multi_treatment_time_corrected, file = "R/object/Model_multi_treatment_time_corrected.RData")

multievent_treatment_corrected2 <- function(){
  
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
    for (tr in c(2,4)){
      phi1_aux[tr,t] <- phi1[tr,1]
      phi2_aux[tr,t] <- phi2[tr,1]
      phi3_aux[tr,t] <- phi3[tr,1]
      
      psi12_aux[tr,t] <- psi12[tr,1]
      psi23_aux[tr,t] <- psi23[tr,1]
      psi13_aux[tr,t] <- psi13[tr,1]
    }
    for (tr in c(1,3)){
      phi1_aux[tr,t] <- phi1[tr+1,1]
      phi2_aux[tr,t] <- phi2[tr+1,1]
      phi3_aux[tr,t] <- phi3[tr+1,1]
      
      psi12_aux[tr,t] <- psi12[tr+1,1]
      psi23_aux[tr,t] <- psi23[tr+1,1]
      psi13_aux[tr,t] <- psi13[tr+1,1]
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

Model_multi_treatment_corrected2 <- jags.parallel(data = jags.data,
                                            inits = inits,
                                            parameters.to.save = parameters,
                                            model.file = multievent_treatment_corrected2,
                                            n.chains = 2,
                                            n.iter = ni)
save(Model_multi_treatment_corrected2, file = "R/object/Model_multi_treatment_corrected2.RData" )


multievent_treatment_corrected3 <- function(){
  
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

Model_multi_treatment_corrected3 <- jags.parallel(data = jags.data,
                                                  inits = inits,
                                                  parameters.to.save = parameters,
                                                  model.file = multievent_treatment_corrected3,
                                                  n.chains = 2,
                                                  n.iter = ni)
save(Model_multi_treatment_corrected3, file = "R/object/Model_multi_treatment_corrected3.RData" )

runtimetot <- Sys.time() - old

print(runtime1,
      runtime2,
      runtime3,
      runtimetot)


#save(Model_multi, file = "R/object/Model_multi.RData")
#save(Model_multi_treatment, file = "R/object/Model_multi_treatment.RData" )
#save(Model_multi_treatment_time, file = "R/object/Model_multi_treatment_time.RData")
